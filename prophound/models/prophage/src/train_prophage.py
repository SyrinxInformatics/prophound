# SPDX-License-Identifier: Apache-2.0
#
# Copyright 2025 Jacob Alford <jalford0000@gmail.com>
#                and Benjamin Alford <ssolar5.dev@gmail.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from __future__ import annotations

import argparse
import json
import os
import random
import signal
import time
from collections import OrderedDict
from dataclasses import asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import yaml
from pyfaidx import Fasta
from sklearn.metrics import average_precision_score, roc_auc_score
from torch.utils.data import DataLoader, Dataset
from tqdm import tqdm

from .config import (  # noqa: F401 -- re-exported for scripts
    DataConfig,
    EarlyStoppingConfig,
    LossConfig,
    ModelConfig,
    ProphageTrainConfig,
    TrainingConfig,
    load_cfg,
)
from .devices import device_supports_amp, resolve_device
from .encoding import collate_batch, onehot_encode, stable_sigmoid
from .model_rc_cnn import build_model

# Utilities

def normalize_contig_id(raw: str) -> List[str]:
    """Return candidate contig IDs for matching manifest values to FASTA headers."""
    out: List[str] = []

    token = raw.split()[0]
    if token:
        out.append(token)

    if "|" in token:
        last = token.split("|")[-1]
        if last and last not in out:
            out.append(last)

    if "." in token:
        base = token.rsplit(".", 1)[0]
        if base and base not in out:
            out.append(base)

    return out


def read_split_ids(splits_dir: str, split: str) -> List[str]:
    path = Path(splits_dir) / f"{split}.txt"
    with open(path, "r") as f:
        return [line.strip() for line in f if line.strip()]


def seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def load_init_weights(
    model: nn.Module,
    checkpoint_path: Path,
    device: torch.device,
    init_mode: str,
) -> None:
    ckpt = torch.load(checkpoint_path, map_location=device, weights_only=False)
    checkpoint_state = ckpt["model"] if isinstance(ckpt, dict) and "model" in ckpt else ckpt

    if init_mode == "strict":
        model.load_state_dict(checkpoint_state, strict=True)
        print(f"[init] loaded model weights from {checkpoint_path} (mode=strict)")
        return

    model_state = model.state_dict()

    if init_mode == "compatible":
        candidate_items = checkpoint_state.items()
    elif init_mode == "encoder_only":
        candidate_items = (
            (name, value)
            for name, value in checkpoint_state.items()
            if name.startswith("encoder.")
        )
    else:
        raise ValueError(
            "training.init_mode must be one of: strict, compatible, encoder_only"
        )

    filtered_state = {
        name: value
        for name, value in candidate_items
        if name in model_state and model_state[name].shape == value.shape
    }
    missing_keys, unexpected_keys = model.load_state_dict(filtered_state, strict=False)
    loaded_count = len(filtered_state)
    print(
        f"[init] loaded {loaded_count} tensors from {checkpoint_path} "
        f"(mode={init_mode}, missing={len(missing_keys)}, unexpected={len(unexpected_keys)})"
    )


def parse_slurm_timelimit(raw: str) -> Optional[int]:
    """Parse `SLURM_TIMELIMIT` into seconds."""
    try:
        raw = raw.strip()
        if "-" in raw:
            days_s, hms = raw.split("-", 1)
            days = int(days_s)
        else:
            days = 0
            hms = raw
        parts = [int(p) for p in hms.split(":")]
        if len(parts) == 3:
            h, m, s = parts
        elif len(parts) == 2:
            h, m, s = 0, parts[0], parts[1]
        else:
            return None
        return days * 86400 + h * 3600 + m * 60 + s
    except Exception:
        return None


# Dataset

class WindowDataset(Dataset):
    """Row-backed window dataset built from the filtered manifest."""

    def __init__(
        self,
        df: pd.DataFrame,
        window_len: int,
    ):
        self.df = df.reset_index(drop=True)
        self.window_len = window_len
        self._fa: "OrderedDict[str, Fasta]" = OrderedDict()
        self._fa_max = 32
        self._fa_index: Dict[str, Dict[str, str]] = {}
        self._missing_contigs: Dict[str, int] = {}

    def __len__(self) -> int:
        return len(self.df)

    def close(self) -> None:
        for _, fa in list(self._fa.items()):
            try:
                fa.close()
            except Exception:
                pass
        self._fa.clear()
        self._fa_index.clear()

    def __del__(self):
        self.close()

    def _get_fa(self, fasta_path: str) -> Fasta:
        fa = self._fa.get(fasta_path)
        if fa is not None:
            self._fa.move_to_end(fasta_path)
            return fa

        fa = Fasta(fasta_path, as_raw=True, sequence_always_upper=True)

        self._fa[fasta_path] = fa
        self._fa.move_to_end(fasta_path)

        # Build a contig index once per FASTA for defensive header matching.
        try:
            keys = list(fa.keys())
        except Exception:
            keys = []
        idx: Dict[str, str] = {}
        for k in keys:
            for cand in normalize_contig_id(k):
                idx.setdefault(cand, k)
        self._fa_index[fasta_path] = idx

        # Keep the FASTA cache bounded so workers do not exhaust file handles.
        while len(self._fa) > self._fa_max:
            _, old = self._fa.popitem(last=False)
            try:
                old.close()
            except Exception:
                pass

        return fa

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        row = self.df.iloc[idx]
        fasta_path = row["fasta_path"]
        contig_id = row["contig_id"]
        start = int(row["start"])
        end = int(row["end"])
        y = float(row["label"])

        fa = self._get_fa(fasta_path)

        # Resolve contig IDs defensively because FASTA headers can vary.
        try_id = contig_id
        if try_id not in fa:
            idx_map = self._fa_index.get(fasta_path, {})
            for cand in normalize_contig_id(contig_id):
                if cand in idx_map:
                    try_id = idx_map[cand]
                    break

        if try_id not in fa:
            key = f"{fasta_path}::{contig_id}"
            self._missing_contigs[key] = self._missing_contigs.get(key, 0) + 1
            if self._missing_contigs[key] <= 3:
                print(
                    f"[warn] contig_id not found in FASTA: {contig_id} ({fasta_path})",
                    flush=True,
                )
            raise KeyError(f"{contig_id} not in {fasta_path}")

        # pyfaidx slices use 0-based half-open coordinates.
        seq = fa[try_id][start:end]
        # Pad or trim only as a safeguard.
        if len(seq) < self.window_len:
            seq = seq + ("N" * (self.window_len - len(seq)))
        elif len(seq) > self.window_len:
            seq = seq[: self.window_len]

        x = onehot_encode(seq)  # (5, L)
        return x, torch.tensor(y, dtype=torch.float32)


# Metrics

@torch.no_grad()
def evaluate(
    model: nn.Module, loader: DataLoader, device: torch.device
) -> Dict[str, float]:
    model.eval()
    all_logits: List[np.ndarray] = []
    all_y: List[np.ndarray] = []

    for batch in tqdm(loader, desc="eval", leave=False):
        x, y = batch
        x = x.to(device, non_blocking=True)
        y = y.to(device, non_blocking=True)
        logits = model(x).view(-1)
        y = y.view(-1)
        all_logits.append(logits.detach().cpu().numpy())
        all_y.append(y.detach().cpu().numpy())

    logits = (
        np.concatenate(all_logits) if all_logits else np.array([], dtype=np.float32)
    )
    y = np.concatenate(all_y) if all_y else np.array([], dtype=np.float32)

    if len(y) == 0:
        return {"auroc": float("nan"), "auprc": float("nan")}

    probs = stable_sigmoid(logits)
    # Some splits can be single-class.
    out = {}
    try:
        out["auroc"] = float(roc_auc_score(y, probs))
    except Exception:
        out["auroc"] = float("nan")
    try:
        out["auprc"] = float(average_precision_score(y, probs))
    except Exception:
        out["auprc"] = float("nan")
    return out


@torch.no_grad()
def predict_logits(
    model: nn.Module, loader: DataLoader, device: torch.device
) -> np.ndarray:
    model.eval()
    all_logits: List[np.ndarray] = []
    for batch in tqdm(loader, desc="predict", leave=False):
        x = batch[0]
        x = x.to(device, non_blocking=True)
        logits = model(x).view(-1)
        all_logits.append(logits.detach().cpu().numpy())
    if not all_logits:
        return np.array([], dtype=np.float32)
    return np.concatenate(all_logits)


# Training

def compute_pos_weight(train_df: pd.DataFrame) -> float:
    n_pos = float((train_df["label"] == 1).sum())
    n_neg = float((train_df["label"] == 0).sum())
    if n_pos <= 0:
        return 1.0
    return n_neg / n_pos


def make_dfs(
    manifest_csv: Optional[str], splits_dir: Optional[str], data_config: DataConfig
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if not manifest_csv or not splits_dir:
        raise ValueError("training and evaluation need data.manifest_csv and data.splits_dir in the config")
    print(f"[data] reading manifest: {manifest_csv}")
    required_columns = [
        "genome_id",
        "contig_id",
        "start",
        "end",
        "label",
        "neg_type",
        "fasta_path",
    ]
    manifest_df = pd.read_csv(manifest_csv, usecols=required_columns)

    manifest_df["label"] = manifest_df["label"].astype(np.int8)
    manifest_df["neg_type"] = manifest_df["neg_type"].astype("category")
    manifest_df["genome_id"] = manifest_df["genome_id"].astype("category")

    split_ids = {
        "train": set(read_split_ids(splits_dir, "train")),
        "val": set(read_split_ids(splits_dir, "val")),
        "test": set(read_split_ids(splits_dir, "test")),
    }

    def subset_split(split_name: str, allowed_neg_types: List[str]) -> pd.DataFrame:
        in_split = manifest_df["genome_id"].isin(split_ids[split_name])
        in_allowed_types = manifest_df["neg_type"].isin(allowed_neg_types)
        return manifest_df[in_split & in_allowed_types].copy()

    train_df = subset_split("train", data_config.use_neg_types_in_train)
    val_df = subset_split("val", data_config.use_neg_types_in_val)
    test_df = subset_split("test", data_config.use_neg_types_in_test)

    print(
        f"[data] train rows: {len(train_df):,} (pos={int((train_df.label==1).sum()):,} neg={int((train_df.label==0).sum()):,})"
    )
    print(
        f"[data]   val rows: {len(val_df):,} (pos={int((val_df.label==1).sum()):,} neg={int((val_df.label==0).sum()):,})"
    )
    print(
        f"[data]  test rows: {len(test_df):,} (pos={int((test_df.label==1).sum()):,} neg={int((test_df.label==0).sum()):,})"
    )

    if len(train_df) == 0:
        raise RuntimeError(
            "Train split produced 0 rows. Likely split IDs don't match manifest genome_id values."
        )
    if int((train_df.label == 1).sum()) == 0:
        raise RuntimeError(
            "Train split has 0 positive windows. Check your splits and BED labeling."
        )
    return train_df, val_df, test_df


def append_eval_metrics(
    history_path: Path, split: str, checkpoint_path: str, metrics: Dict[str, float]
) -> None:
    if history_path.exists():
        with open(history_path, "r") as f:
            history = json.load(f)
    else:
        history = {"train": [], "val": [], "test": []}
    history.setdefault(split, [])
    history[split].append(
        {
            "epoch": "eval_only",
            "checkpoint": checkpoint_path,
            **metrics,
        }
    )
    with open(history_path, "w") as f:
        json.dump(history, f, indent=2)


def eval_only(
    cfg: ProphageTrainConfig,
    out_dir: Path,
    checkpoint_path: str,
    split: str,
    requested_device: str,
) -> None:
    device = resolve_device(requested_device)
    print(f"[sys] device: {device} (requested={requested_device})")
    print(f"[eval] split={split} checkpoint={checkpoint_path}")

    train_df, val_df, test_df = make_dfs(
        cfg.data.manifest_csv, cfg.data.splits_dir, cfg.data
    )
    split_df_map = {
        "train": train_df,
        "val": val_df,
        "test": test_df,
    }
    target_df = split_df_map[split]
    if len(target_df) == 0:
        raise RuntimeError(f"{split} split produced 0 rows.")

    target_ds = WindowDataset(target_df, window_len=cfg.data.window_len)
    target_loader = DataLoader(
        target_ds,
        batch_size=cfg.training.batch_size,
        shuffle=False,
        num_workers=cfg.training.num_workers,
        pin_memory=True,
        collate_fn=collate_batch,
        persistent_workers=False,
    )

    model = build_model(cfg, device)
    ckpt = torch.load(checkpoint_path, map_location=device, weights_only=False)
    state_dict = ckpt["model"] if isinstance(ckpt, dict) and "model" in ckpt else ckpt
    model.load_state_dict(state_dict, strict=True)

    metrics = evaluate(model, target_loader, device)
    print(
        f"[eval:{split}] auroc={metrics['auroc']:.4f} auprc={metrics['auprc']:.4f}",
        flush=True,
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "metrics").mkdir(exist_ok=True)
    append_eval_metrics(
        out_dir / "metrics" / "history.json",
        split=split,
        checkpoint_path=checkpoint_path,
        metrics=metrics,
    )
    target_ds.close()


def export_predictions(
    cfg: ProphageTrainConfig,
    out_dir: Path,
    checkpoint_path: str,
    split: str,
    out_csv: str,
    requested_device: str,
) -> None:
    device = resolve_device(requested_device)
    print(f"[sys] device: {device} (requested={requested_device})")
    print(
        f"[predict] split={split} checkpoint={checkpoint_path} out_csv={out_csv}",
        flush=True,
    )

    train_df, val_df, test_df = make_dfs(
        cfg.data.manifest_csv, cfg.data.splits_dir, cfg.data
    )
    split_df_map = {
        "train": train_df,
        "val": val_df,
        "test": test_df,
    }
    target_df = split_df_map[split].reset_index(drop=True)
    if len(target_df) == 0:
        raise RuntimeError(f"{split} split produced 0 rows.")

    target_ds = WindowDataset(target_df, window_len=cfg.data.window_len)
    target_loader = DataLoader(
        target_ds,
        batch_size=cfg.training.batch_size,
        shuffle=False,
        num_workers=cfg.training.num_workers,
        pin_memory=True,
        collate_fn=collate_batch,
        persistent_workers=False,
    )

    model = build_model(cfg, device)
    ckpt = torch.load(checkpoint_path, map_location=device, weights_only=False)
    state_dict = ckpt["model"] if isinstance(ckpt, dict) and "model" in ckpt else ckpt
    model.load_state_dict(state_dict, strict=True)

    logits = predict_logits(model, target_loader, device)
    if len(logits) != len(target_df):
        raise RuntimeError(
            f"Prediction row mismatch: logits={len(logits)} target_df={len(target_df)}"
        )
    probs = stable_sigmoid(logits)
    preds = (probs >= 0.5).astype(np.int8)

    out_df = target_df[
        ["genome_id", "contig_id", "start", "end", "label", "neg_type"]
    ].copy()
    out_df["checkpoint"] = checkpoint_path
    out_df["split"] = split
    out_df["logit"] = logits
    out_df["prob"] = probs
    out_df["pred_label"] = preds

    out_path = Path(out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, index=False)
    print(f"[predict:{split}] wrote {len(out_df):,} rows -> {out_path}", flush=True)
    target_ds.close()


def train_one_run(cfg: ProphageTrainConfig, out_dir: Path, requested_device: str) -> None:
    seed_everything(cfg.training.seed)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "checkpoints").mkdir(exist_ok=True)
    (out_dir / "metrics").mkdir(exist_ok=True)

    device = resolve_device(requested_device)
    print(f"[sys] device: {device} (requested={requested_device})")

    train_df, val_df, test_df = make_dfs(
        cfg.data.manifest_csv, cfg.data.splits_dir, cfg.data
    )

    # pos_weight
    if cfg.loss.pos_weight == "auto":
        pos_weight_value = compute_pos_weight(train_df)
    else:
        pos_weight_value = float(cfg.loss.pos_weight)
    print(f"[loss] pos_weight: {pos_weight_value:.4f}")
    pos_weight = torch.tensor([pos_weight_value], dtype=torch.float32, device=device)

    # datasets/loaders
    train_ds = WindowDataset(train_df, window_len=cfg.data.window_len)
    val_ds = (
        WindowDataset(val_df, window_len=cfg.data.window_len) if len(val_df) else None
    )
    test_ds = (
        WindowDataset(test_df, window_len=cfg.data.window_len) if len(test_df) else None
    )

    train_loader = DataLoader(
        train_ds,
        batch_size=cfg.training.batch_size,
        shuffle=True,
        num_workers=cfg.training.num_workers,
        pin_memory=True,
        collate_fn=collate_batch,
        persistent_workers=False,
    )
    val_loader = (
        DataLoader(
            val_ds,
            batch_size=cfg.training.batch_size,
            shuffle=False,
            num_workers=cfg.training.num_workers,
            pin_memory=True,
            collate_fn=collate_batch,
            persistent_workers=False,
        )
        if val_ds is not None
        else None
    )

    # model
    model = build_model(cfg, device)

    optimizer = torch.optim.AdamW(
        model.parameters(), lr=cfg.training.lr, weight_decay=cfg.training.weight_decay
    )
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    amp_enabled = cfg.training.amp and device_supports_amp(device)
    scaler = torch.amp.GradScaler("cuda", enabled=amp_enabled)
    print(
        "[train] "
        f"lr={cfg.training.lr} grad_clip={cfg.training.grad_clip} amp={amp_enabled} "
        f"workers={cfg.training.num_workers} "
        f"save_every_sec={cfg.training.save_every_sec} "
        f"save_every_steps={cfg.training.save_every_steps}",
        flush=True,
    )

    if cfg.training.early_stopping.mode == "min":
        best_score = float("inf")
    else:
        best_score = -1e9
    best_path = out_dir / "checkpoints" / "best.pt"
    last_path = out_dir / "checkpoints" / "last.pt"

    patience = cfg.training.early_stopping.patience
    mode = cfg.training.early_stopping.mode
    metric_name = cfg.training.early_stopping.metric
    bad_epochs = 0

    history: Dict[str, List[dict]] = {"train": [], "val": [], "test": []}
    nonfinite_loss_skips = 0
    nonfinite_grad_skips = 0

    cfg_dict = asdict(cfg)
    # snapshot config for provenance
    cfg_path = out_dir / "config_snapshot.yaml"
    if not cfg_path.exists():
        with open(cfg_path, "w") as f:
            yaml.safe_dump(cfg_dict, f)

    # resume if available
    start_epoch = 1
    resume_step = 0
    global_step = 0
    bad_epochs = 0

    if cfg.training.auto_resume and last_path.exists():
        ckpt = torch.load(last_path, map_location=device, weights_only=False)
        model.load_state_dict(ckpt["model"])
        if "opt" in ckpt:
            optimizer.load_state_dict(ckpt["opt"])
        if "scaler" in ckpt and ckpt["scaler"] is not None:
            try:
                scaler.load_state_dict(ckpt["scaler"])
            except Exception:
                pass
        ckpt_epoch = int(ckpt.get("epoch", 0))
        ckpt_epoch_step = int(ckpt.get("epoch_step", 0))
        if "epoch_complete" in ckpt:
            if bool(ckpt.get("epoch_complete")):
                start_epoch = ckpt_epoch + 1
                resume_step = 0
            else:
                start_epoch = ckpt_epoch
                resume_step = ckpt_epoch_step
        else:
            # Backward compatibility for older checkpoints: prior behavior.
            start_epoch = ckpt_epoch + 1
            resume_step = 0
        global_step = int(ckpt.get("step", 0))
        best_score = float(ckpt.get("best_score", best_score))
        bad_epochs = int(ckpt.get("bad_epochs", 0))
        history = ckpt.get("history", history)

        # restore RNG
        rng = ckpt.get("rng", {})
        try:
            random.setstate(rng.get("python"))
            np.random.set_state(rng.get("numpy"))
            torch.set_rng_state(rng.get("torch"))
            if device.type == "cuda" and "cuda" in rng:
                torch.cuda.set_rng_state_all(rng.get("cuda"))
        except Exception:
            pass

        print(
            f"[ckpt] resumed from {last_path} at epoch {start_epoch} "
            f"(resume_step={resume_step}) global_step={global_step}"
        )
    elif cfg.training.init_checkpoint:
        init_path = Path(cfg.training.init_checkpoint)
        if init_path.exists():
            load_init_weights(
                model,
                checkpoint_path=init_path,
                device=device,
                init_mode=cfg.training.init_mode,
            )
        else:
            raise FileNotFoundError(
                f"training.init_checkpoint not found: {init_path}"
            )

    # walltime awareness
    job_start = time.time()
    slurm_limit = parse_slurm_timelimit(os.environ.get("SLURM_TIMELIMIT", ""))
    slurm_start_raw = os.environ.get("SLURM_JOB_START_TIME")
    if slurm_start_raw and slurm_start_raw.isdigit():
        job_start = float(slurm_start_raw)

    stop_training = {"flag": False}

    def handle_sigterm(signum, frame):
        stop_training["flag"] = True

    signal.signal(signal.SIGTERM, handle_sigterm)
    signal.signal(signal.SIGINT, handle_sigterm)

    def score_from(metrics: Dict[str, float]) -> float:
        if metric_name.endswith("auprc"):
            return metrics.get("auprc", float("nan"))
        if metric_name.endswith("auroc"):
            return metrics.get("auroc", float("nan"))
        return metrics.get("auprc", float("nan"))

    last_ckpt_time = time.time()

    def save_checkpoint(
        path: Path,
        epoch: int,
        step: int,
        epoch_step: int = 0,
        epoch_complete: bool = False,
    ) -> None:
        rng = {
            "python": random.getstate(),
            "numpy": np.random.get_state(),
            "torch": torch.get_rng_state(),
        }
        if device.type == "cuda":
            rng["cuda"] = torch.cuda.get_rng_state_all()

        torch.save(
            {
                "model": model.state_dict(),
                "opt": optimizer.state_dict(),
                "scaler": scaler.state_dict() if scaler is not None else None,
                "epoch": epoch,
                "step": step,
                "epoch_step": epoch_step,
                "epoch_complete": epoch_complete,
                "best_score": best_score,
                "bad_epochs": bad_epochs,
                "history": history,
                "cfg": cfg_dict,
                "rng": rng,
            },
            path,
        )

    for epoch in range(start_epoch, cfg.training.epochs + 1):
        model.train()
        epoch_start_time = time.time()
        running_loss = 0.0
        seen_examples = 0
        epoch_resume_step = resume_step if epoch == start_epoch else 0

        pbar = tqdm(
            train_loader, desc=f"epoch {epoch}/{cfg.training.epochs}", leave=True
        )
        for step, batch in enumerate(pbar, start=1):
            if epoch_resume_step and step <= epoch_resume_step:
                continue

            x, y = batch
            x = x.to(device, non_blocking=True)
            y = y.to(device, non_blocking=True).view(-1)

            optimizer.zero_grad(set_to_none=True)
            with torch.amp.autocast("cuda", enabled=amp_enabled):
                logits = model(x).view(-1)
                loss = criterion(logits, y)

            if not torch.isfinite(loss):
                nonfinite_loss_skips += 1
                print(
                    f"[warn] non-finite loss epoch={epoch} step={step} loss={loss.item()} "
                    f"(skipped={nonfinite_loss_skips})",
                    flush=True,
                )
                optimizer.zero_grad(set_to_none=True)
                continue

            if scaler is not None and scaler.is_enabled():
                scaler.scale(loss).backward()
            else:
                loss.backward()

            if step == 1 or step % 1000 == 0:
                print(
                    f"[epoch {epoch}] step {step}/{len(train_loader)} loss={loss.item():.5f}",
                    flush=True,
                )

            if cfg.training.grad_clip and cfg.training.grad_clip > 0:
                if scaler is not None and scaler.is_enabled():
                    scaler.unscale_(optimizer)
                nn.utils.clip_grad_norm_(model.parameters(), cfg.training.grad_clip)

            grads_finite = True
            for p in model.parameters():
                if p.grad is None:
                    continue
                if not torch.isfinite(p.grad).all():
                    grads_finite = False
                    break
            if not grads_finite:
                nonfinite_grad_skips += 1
                print(
                    f"[warn] non-finite gradients epoch={epoch} step={step} "
                    f"(skipped={nonfinite_grad_skips})",
                    flush=True,
                )
                optimizer.zero_grad(set_to_none=True)
                if scaler is not None and scaler.is_enabled():
                    scaler.update()
                continue

            if scaler is not None and scaler.is_enabled():
                scaler.step(optimizer)
                scaler.update()
            else:
                optimizer.step()

            global_step += 1

            if cfg.training.save_every_sec > 0 and (time.time() - last_ckpt_time) >= cfg.training.save_every_sec:
                save_checkpoint(
                    last_path,
                    epoch,
                    global_step,
                    epoch_step=step,
                    epoch_complete=False,
                )
                last_ckpt_time = time.time()
                print(f"[ckpt] wrote last.pt at epoch {epoch} step {global_step}", flush=True)
            elif cfg.training.save_every_steps > 0 and (global_step % cfg.training.save_every_steps) == 0:
                save_checkpoint(
                    last_path,
                    epoch,
                    global_step,
                    epoch_step=step,
                    epoch_complete=False,
                )
                print(
                    f"[ckpt] wrote last.pt at epoch {epoch} step {global_step} "
                    f"(step interval)",
                    flush=True,
                )

            # walltime cutoff
            if slurm_limit is not None:
                elapsed = time.time() - job_start
                remaining = slurm_limit - elapsed
                if remaining <= cfg.training.slurm_margin_sec:
                    save_checkpoint(
                        last_path,
                        epoch,
                        global_step,
                        epoch_step=step,
                        epoch_complete=False,
                    )
                    print("[ckpt] nearing walltime; saved and exiting.", flush=True)
                    train_ds.close()
                    if val_ds is not None:
                        val_ds.close()
                    if test_ds is not None:
                        test_ds.close()
                    return

            if stop_training["flag"]:
                save_checkpoint(
                    last_path,
                    epoch,
                    global_step,
                    epoch_step=step,
                    epoch_complete=False,
                )
                print("[ckpt] received termination signal; saved and exiting.", flush=True)
                train_ds.close()
                if val_ds is not None:
                    val_ds.close()
                if test_ds is not None:
                    test_ds.close()
                return

            running_loss += float(loss.item()) * y.size(0)
            seen_examples += y.size(0)
            pbar.set_postfix(loss=(running_loss / max(seen_examples, 1)))

        train_loss = running_loss / max(seen_examples, 1)
        epoch_time = time.time() - epoch_start_time

        train_metrics = {"loss": float(train_loss), "time_sec": float(epoch_time)}
        history["train"].append({"epoch": epoch, **train_metrics})

        print(f"[epoch {epoch}] train loss={train_loss:.5f} time={epoch_time:.1f}s")
        print(
            f"[epoch {epoch}] skipped batches: nonfinite_loss={nonfinite_loss_skips} "
            f"nonfinite_grad={nonfinite_grad_skips}",
            flush=True,
        )

        # validation
        if val_loader is not None and val_ds is not None and len(val_ds) > 0:
            val_metrics = evaluate(model, val_loader, device)
            history["val"].append({"epoch": epoch, **val_metrics})
            print(
                f"[epoch {epoch}] val auroc={val_metrics['auroc']:.4f} auprc={val_metrics['auprc']:.4f}"
            )
            validation_score = score_from(val_metrics)
        else:
            validation_score = float("nan")

        # checkpoint
        save_checkpoint(last_path, epoch, global_step, epoch_step=0, epoch_complete=True)

        improved = False
        if not np.isnan(validation_score):
            if mode == "max" and validation_score > best_score:
                improved = True
            elif mode == "min" and validation_score < best_score:
                improved = True

        if improved:
            best_score = validation_score
            bad_epochs = 0
            save_checkpoint(best_path, epoch, global_step, epoch_step=0, epoch_complete=True)
            print(f"[ckpt] new best {metric_name}={best_score:.4f} -> {best_path}")
        else:
            bad_epochs += 1
            print(f"[early] bad_epochs={bad_epochs}/{patience}")
            if bad_epochs >= patience:
                print("[early] stopping")
                break

        # write metrics each epoch
        with open(out_dir / "metrics" / "history.json", "w") as f:
            json.dump(history, f, indent=2)

    # final test eval (best)
    if best_path.exists() and test_ds is not None and len(test_ds) > 0:
        ckpt = torch.load(best_path, map_location=device, weights_only=False)
        model.load_state_dict(ckpt["model"])
        test_loader = DataLoader(
            test_ds,
            batch_size=cfg.training.batch_size,
            shuffle=False,
            num_workers=cfg.training.num_workers,
            pin_memory=True,
            collate_fn=collate_batch,
            persistent_workers=False,
        )
        test_metrics = evaluate(model, test_loader, device)
        history["test"].append({"epoch": "best", **test_metrics})
        print(
            f"[test] auroc={test_metrics['auroc']:.4f} auprc={test_metrics['auprc']:.4f}"
        )
        with open(out_dir / "metrics" / "history.json", "w") as f:
            json.dump(history, f, indent=2)

    train_ds.close()
    if val_ds is not None:
        val_ds.close()
    if test_ds is not None:
        test_ds.close()
    print("[done]")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--config", required=True, help="Path to runs/prophage_v0/config.yaml"
    )
    ap.add_argument(
        "--eval-only",
        action="store_true",
        help="Run evaluation only on a saved checkpoint.",
    )
    ap.add_argument(
        "--checkpoint",
        default=None,
        help="Checkpoint path for --eval-only. Defaults to runs/<run_name>/checkpoints/best.pt.",
    )
    ap.add_argument(
        "--split",
        choices=["train", "val", "test"],
        default="test",
        help="Split to evaluate for --eval-only.",
    )
    ap.add_argument(
        "--export-predictions",
        action="store_true",
        help="Export per-window predictions for a saved checkpoint.",
    )
    ap.add_argument(
        "--predictions-out",
        default=None,
        help="CSV path for --export-predictions. Defaults to runs/<run_name>/metrics/<split>_predictions.csv.",
    )
    ap.add_argument(
        "--device",
        default="auto",
        choices=["auto", "cpu", "cuda", "rocm", "mps", "xpu", "intel"],
        help="Inference/training device. 'auto' prefers CUDA, then MPS, then XPU, then CPU.",
    )
    args = ap.parse_args()

    cfg = load_cfg(args.config)
    out_dir = Path("runs") / cfg.run_name
    print(f"[run] {cfg.run_name}")
    print(f"[out] {out_dir}")
    if args.eval_only and args.export_predictions:
        raise ValueError("Use either --eval-only or --export-predictions, not both.")
    if args.eval_only:
        checkpoint_path = args.checkpoint or str(out_dir / "checkpoints" / "best.pt")
        eval_only(
            cfg,
            out_dir,
            checkpoint_path=checkpoint_path,
            split=args.split,
            requested_device=args.device,
        )
    elif args.export_predictions:
        checkpoint_path = args.checkpoint or str(out_dir / "checkpoints" / "best.pt")
        predictions_out = args.predictions_out or str(
            out_dir / "metrics" / f"{args.split}_predictions.csv"
        )
        export_predictions(
            cfg,
            out_dir,
            checkpoint_path=checkpoint_path,
            split=args.split,
            out_csv=predictions_out,
            requested_device=args.device,
        )
    else:
        train_one_run(cfg, out_dir, requested_device=args.device)


if __name__ == "__main__":
    main()
