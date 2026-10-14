#!/usr/bin/env python3
"""Sanity checks for SurrealGAN training.

SurrealGAN's patient batch iterator runs int(n_patients / batchsize) - 1 training steps per
epoch. When n_patients < 2 * batchsize that is zero steps, so the networks never train: every
saved epoch holds the same random weights and `Rindices_corr` in model_agreements.csv is
identical at every epoch. Depending on the installed version, such a run either crashes at
epoch 100 or saves the untrained networks without any warning.

Two uses:

  # before submitting: training steps per epoch, and the largest usable batch size
  python check_model_training.py --data dt_Model_learning_1014.csv --batchsize 300

  # after training: scan output folders (a configuration folder, or a parent holding Model_*)
  python check_model_training.py /path/to/Multidisease

The scan reports, per configuration folder, how much `Rindices_corr` varies across the saved
epochs and whether the saved weights differ between the first and last epoch of one model.
A folder is FROZEN when the agreement does not change, or the weights are identical. The exit
status is 1 if any folder is FROZEN. The weight check needs PyTorch; without it only the CSV
is checked. The steps check needs SurrealGAN, because it counts batches with the installed
package's own iterator.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ID_COLUMNS = ["participant_id", "diagnosis"]
# Rindices_corr differs by about 1e-2 between epochs of a training model, and by about 1e-17
# (floating-point noise) between epochs of a model that is not training.
FROZEN_TOLERANCE = 1e-9
MIN_ROWS = 3
NETWORKS = ("netDecomposer", "netReconstruction")


def patient_matrix(data):
    """Feature matrix of the patients (diagnosis == 1), as SurrealGAN's training iterator gets it."""
    return data.loc[data["diagnosis"] == 1].drop(columns=ID_COLUMNS).to_numpy()


def steps_per_epoch(data, batchsize, fraction=1.0):
    """Training steps SurrealGAN runs per epoch, counted with the installed PTIterator."""
    from SurrealGAN.data_loading import PTIterator

    state = np.random.get_state()          # PTIterator seeds numpy's global generator
    try:
        return sum(1 for _ in PTIterator(patient_matrix(data), 0, fraction, batchsize))
    finally:
        np.random.set_state(state)


def largest_batchsize(data, min_steps, fraction=1.0):
    """Largest batch size that gives at least min_steps steps per epoch, or None."""
    for batchsize in range(int(fraction * len(patient_matrix(data))), 0, -1):
        if steps_per_epoch(data, batchsize, fraction) >= min_steps:
            return batchsize
    return None


def agreement_summary(csv_path):
    """(rows, distinct epochs, max - min of Rindices_corr) of a model_agreements.csv."""
    table = pd.read_csv(csv_path)
    values = pd.to_numeric(table["Rindices_corr"], errors="coerce")
    spread = float(values.max() - values.min()) if len(values) else float("nan")
    return len(table), int(table["epoch"].nunique()), spread


def weights_change(model_file):
    """True if the weights differ between the first and last saved epoch, False if identical,
    None if fewer than two epochs are saved or PyTorch is missing."""
    try:
        import torch
    except ImportError:
        return None
    try:
        checkpoint = torch.load(model_file, weights_only=False)
    except TypeError:                      # PyTorch before 1.13 has no weights_only argument
        checkpoint = torch.load(model_file)
    epochs = sorted(key for key in checkpoint if isinstance(key, int))
    if len(epochs) < 2:
        return None
    first, last = checkpoint[epochs[0]], checkpoint[epochs[-1]]
    return any(not torch.equal(tensor, last[net][name])
               for net in NETWORKS for name, tensor in first[net].items())


def configuration_folders(path):
    """The folder itself if it holds model_agreements.csv, else its Model_* subfolders."""
    if (path / "model_agreements.csv").is_file():
        return [path]
    return sorted(p for p in path.glob("Model_*") if p.is_dir())


def check_folder(folder):
    """(rows, spread, weights, verdict) for one configuration folder."""
    csv_path = folder / "model_agreements.csv"
    if not csv_path.is_file():
        return None, None, None, "NO model_agreements.csv"
    rows, epochs, spread = agreement_summary(csv_path)
    model_files = sorted(folder.glob("model[0-9]*"))
    weights = weights_change(model_files[0]) if model_files else None
    if weights is False or (rows >= MIN_ROWS and spread < FROZEN_TOLERANCE):
        verdict = "FROZEN"
    elif epochs < MIN_ROWS and weights is None:
        verdict = "too few epochs to judge"
    else:
        verdict = "ok"
    return rows, spread, weights, verdict


def main():
    parser = argparse.ArgumentParser(description="Sanity checks for SurrealGAN training")
    parser.add_argument("folders", nargs="*", type=Path,
                        help="output folders: a configuration folder, or a parent holding Model_* folders")
    parser.add_argument("--data", type=Path, help="training CSV (participant_id, diagnosis, features)")
    parser.add_argument("--batchsize", type=int, help="batch size to check with --data")
    parser.add_argument("--fraction", type=float, default=1.0, help="data fraction per repetition (default 1)")
    args = parser.parse_args()
    if not args.folders and not args.data:
        parser.error("give output folders to scan, or --data and --batchsize")
    failed = False

    if args.data:
        data = pd.read_csv(args.data)
        n_patients = int((data["diagnosis"] == 1).sum())
        print(f"{n_patients} patients in {args.data.name}")
        if args.batchsize:
            steps = steps_per_epoch(data, args.batchsize, args.fraction)
            print(f"batchsize {args.batchsize}: {steps} training steps per epoch")
            failed |= steps < 1
        for minimum in (1, 3):
            print(f"largest batch size with at least {minimum} step(s) per epoch: "
                  f"{largest_batchsize(data, minimum, args.fraction)}")
        if failed:
            print("FROZEN: this batch size gives no training steps, so the model would never train")

    for root in args.folders:
        folders = configuration_folders(root)
        if not folders:
            print(f"{root}: no configuration folders found")
            continue
        print(f"{'folder':<34} {'rows':>5} {'Rindices_corr range':>20} {'weights change':>15}  verdict")
        for folder in folders:
            rows, spread, weights, verdict = check_folder(folder)
            failed |= verdict == "FROZEN"
            shown = "n/a" if weights is None else str(weights)
            spread_text = "n/a" if spread is None else f"{spread:.3g}"
            print(f"{folder.name:<34} {'' if rows is None else rows:>5} {spread_text:>20} {shown:>15}  {verdict}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
