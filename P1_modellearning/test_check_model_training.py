"""Tests of check_model_training.py (needs pandas; the step and weight checks need SurrealGAN / PyTorch).

The step counts below describe the released SurrealGAN (0.1.1 and current main), whose patient
iterator runs int(n_patients / batchsize) - 1 steps per epoch. If an upstream release changes
that, these expectations change with it; the code under test counts with the installed package.
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import check_model_training as cmt

HERE = Path(__file__).resolve().parent


def make_data(n_patients=411, n_controls=421, n_features=5):
    rng = np.random.default_rng(0)
    n = n_patients + n_controls
    data = pd.DataFrame(rng.normal(size=(n, n_features)), columns=[f"f{i}" for i in range(n_features)])
    data.insert(0, "diagnosis", [-1] * n_controls + [1] * n_patients)
    data.insert(0, "participant_id", [f"S{i}" for i in range(n)])
    return data


def write_agreements(folder, values):
    folder.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"epoch": [4000 * (i + 1) for i in range(len(values))], "Rindices_corr": values}).to_csv(
        folder / "model_agreements.csv", index=False)


def test_steps_per_epoch_counts_with_the_installed_iterator():
    pytest.importorskip("SurrealGAN")
    data = make_data()
    assert cmt.steps_per_epoch(data, 300) == 0          # the situation in P1_Model_test.py: no training
    assert cmt.steps_per_epoch(data, 200) == 1
    assert cmt.steps_per_epoch(data, 100) == 3
    assert cmt.steps_per_epoch(data, 411) == 0


def test_largest_batchsize_and_random_state_untouched():
    pytest.importorskip("SurrealGAN")
    data = make_data()
    state = np.random.get_state()[1].copy()
    assert cmt.largest_batchsize(data, 1) == 205         # int(411/205) = 2 -> 1 step; 206 -> 0 steps
    assert cmt.largest_batchsize(data, 3) == 102         # int(411/102) = 4 -> 3 steps; 103 -> 2 steps
    assert np.array_equal(np.random.get_state()[1], state)
    assert cmt.largest_batchsize(data, 10 ** 6) is None


def test_constant_agreement_is_frozen_and_changing_agreement_is_not(tmp_path):
    write_agreements(tmp_path / "Model_frozen", [0.4360770713885283] * 33 + [0.43607707138852836])
    write_agreements(tmp_path / "Model_training", [0.627, 0.710, 0.765, 0.831, 0.834])
    write_agreements(tmp_path / "Model_short", [0.5, 0.5])
    assert cmt.check_folder(tmp_path / "Model_frozen")[3] == "FROZEN"
    assert cmt.check_folder(tmp_path / "Model_training")[3] == "ok"
    assert cmt.check_folder(tmp_path / "Model_short")[3] == "too few epochs to judge"
    assert cmt.check_folder(tmp_path / "Model_absent")[3] == "NO model_agreements.csv"
    assert [p.name for p in cmt.configuration_folders(tmp_path)] == ["Model_frozen", "Model_short", "Model_training"]
    assert [p.name for p in cmt.configuration_folders(tmp_path / "Model_frozen")] == ["Model_frozen"]


def checkpoint(decomposer_values, reconstruction_values):
    torch = pytest.importorskip("torch")
    return {"nROI": 5,
            4000: {"netDecomposer": {"w": torch.tensor(decomposer_values[0])},
                   "netReconstruction": {"w": torch.tensor(reconstruction_values[0])}},
            8000: {"netDecomposer": {"w": torch.tensor(decomposer_values[1])},
                   "netReconstruction": {"w": torch.tensor(reconstruction_values[1])}}}


def test_weight_check_detects_identical_and_changing_weights(tmp_path):
    torch = pytest.importorskip("torch")
    same = checkpoint(([1.0, 2.0], [1.0, 2.0]), ([0.5], [0.5]))
    changed = checkpoint(([1.0, 2.0], [1.0, 2.5]), ([0.5], [0.5]))
    one_epoch = {"nROI": 5, 4000: same[4000]}
    for name, ck in (("same", same), ("changed", changed), ("one", one_epoch)):
        torch.save(ck, tmp_path / name)
    assert cmt.weights_change(tmp_path / "same") is False
    assert cmt.weights_change(tmp_path / "changed") is True
    assert cmt.weights_change(tmp_path / "one") is None


def test_identical_weights_are_frozen_even_when_the_csv_varies(tmp_path):
    torch = pytest.importorskip("torch")
    folder = tmp_path / "Model_M2_Lam0.1_Gam0.1"
    write_agreements(folder, [0.5, 0.6, 0.7, 0.8])
    torch.save(checkpoint(([1.0, 2.0], [1.0, 2.0]), ([0.5], [0.5])), folder / "model0")
    rows, spread, weights, verdict = cmt.check_folder(folder)
    assert (rows, weights, verdict) == (4, False, "FROZEN") and spread > 0.2


def test_command_line_exit_status(tmp_path):
    pytest.importorskip("SurrealGAN")
    data_path = tmp_path / "data.csv"
    make_data().to_csv(data_path, index=False)
    script = str(HERE / "check_model_training.py")
    bad = subprocess.run([sys.executable, script, "--data", str(data_path), "--batchsize", "300"],
                         capture_output=True, text=True)
    assert bad.returncode == 1 and "0 training steps per epoch" in bad.stdout
    assert "largest batch size with at least 3 step(s) per epoch: 102" in bad.stdout
    good = subprocess.run([sys.executable, script, "--data", str(data_path), "--batchsize", "100"],
                          capture_output=True, text=True)
    assert good.returncode == 0 and "3 training steps per epoch" in good.stdout
    write_agreements(tmp_path / "Model_frozen", [0.4] * 5)
    scan = subprocess.run([sys.executable, script, str(tmp_path)], capture_output=True, text=True)
    assert scan.returncode == 1 and "FROZEN" in scan.stdout
