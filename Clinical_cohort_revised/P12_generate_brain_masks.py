#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_WORKLIST = HERE / "worklist.csv"
DEFAULT_SAMPLE = Path("/cwStorage/nodecw_group/FY_data/QSM_HUASHAN/sample")
REQUIRED_COLUMNS = ("image_dir_id", "hospital_id", "qsm_file", "include_main_model", "runnable")
# Handles both ..._M_e1 and ..._M_P180924_e1 / ..._M_P162013_e3.4.
# The numeric suffix is either an echo index or, in some exports, EchoTime in ms;
# within one series either convention sorts the echoes from earliest to latest.
MAG_ECHO_RE = re.compile(
    r"(?:^|_)M(?:_P\d+[A-Za-z]?(?:-\d+)?)?_e(?P<echo>\d+(?:\.\d+)?)$",
    re.IGNORECASE,
)


def nii_stem(path: Path) -> str:
    name = path.name
    if name.lower().endswith(".nii.gz"):
        return name[:-7]
    if name.lower().endswith(".nii"):
        return name[:-4]
    return path.stem


def sidecar_echo_time(image_path: Path):
    """Read JSON EchoTime when present; return None when absent/unparseable."""
    sidecar = image_path.with_name(nii_stem(image_path) + ".json")
    if not sidecar.is_file():
        return None
    try:
        with sidecar.open("r", encoding="utf-8-sig") as handle:
            metadata = json.load(handle)
        value = metadata.get("EchoTime")
        if isinstance(value, (int, float)):
            return float(value)
    except (OSError, ValueError, TypeError):
        pass
    return None


def read_worklist(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = set(REQUIRED_COLUMNS).difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Worklist is missing columns: {', '.join(sorted(missing))}")
        rows = [row for row in reader
                if row.get("include_main_model", "").strip().lower() == "yes"
                and row.get("runnable", "").strip().lower() == "yes"]

    grouped = defaultdict(list)
    for row in rows:
        iid = row.get("image_dir_id", "").strip()
        if iid:
            grouped[iid].append(row)
    return grouped


def choose_source(iid: str, rows):
    """Find same-directory magnitude echoes and choose the earliest echo."""
    participants = {r.get("hospital_id", "").strip() for r in rows
                    if r.get("hospital_id", "").strip()}
    qsm_paths = {Path(r.get("qsm_file", "").strip()) for r in rows
                 if r.get("qsm_file", "").strip()}
    if len(participants) > 1:
        return None, "", 0, "", f"Image directory ID maps to multiple hospital IDs: {sorted(participants)}"
    if len(qsm_paths) != 1:
        return None, "", 0, "", f"QSM source path is not unique: {len(qsm_paths)} paths"
    qsm = next(iter(qsm_paths))
    if not qsm.is_file():
        return None, "", 0, "", f"QSM source file is missing or unreadable: {qsm}"

    try:
        magnitudes = []
        for path in qsm.parent.iterdir():
            if not path.is_file() or not path.name.lower().endswith((".nii", ".nii.gz")):
                continue
            stem = nii_stem(path)
            match = MAG_ECHO_RE.search(stem)
            if not match or "qsm" not in stem.lower():
                continue
            low = stem.lower()
            if any(token in low for token in ("_ph", "phase", "localp", "deltap", "swi", "swan")):
                continue
            suffix_echo = float(match.group("echo"))
            json_echo = sidecar_echo_time(path)
            magnitudes.append((path, suffix_echo, json_echo))
    except OSError as exc:
        return None, "", 0, "", f"Cannot read the MR directory containing the QSM file {qsm.parent}: {exc}"
    if not magnitudes:
        return None, "", 0, "", "No multi-echo QSM magnitude images matching the naming rule were found in the same directory"

    # Prefer JSON EchoTime only when it is available for every echo; otherwise
    # sort by the numeric suffix (e1/e2 or TE-like e3.4/e6.4) consistently.
    use_json_te = all(item[2] is not None for item in magnitudes)
    key_index = 2 if use_json_te else 1
    ordered = sorted(magnitudes, key=lambda item: (item[key_index], item[0].name.lower()))
    basis = "JSON EchoTime" if use_json_te else "filename echo suffix"
    first_value = ordered[0][key_index]
    first_echoes = [item for item in ordered if item[key_index] == first_value]
    if len(first_echoes) != 1:
        return (None, ";".join(x[0].name for x in ordered), len(ordered),
                basis,
                f"Multiple earliest-echo candidates; cannot select one uniquely: {[x[0].name for x in first_echoes]}")
    selected = ordered[0]
    return selected[0], ";".join(x[0].name for x in ordered), len(ordered), basis, ""


def run_bet(iid: str, magnitude: Path, sample_root: Path,
            fraction: float, overwrite: bool):
    out_dir = sample_root / iid
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / f"{iid}_M_e1_brain"
    mask = out_dir / f"{iid}_M_e1_brain_mask.nii.gz"
    brain = out_dir / f"{iid}_M_e1_brain.nii.gz"

    if mask.is_file() and mask.stat().st_size > 0 and not overwrite:
        return "EXISTS", str(mask), "Existing mask skipped (use --overwrite to regenerate)"
    if not overwrite and brain.exists():
        return "CONFLICT", str(mask), f"Brain-extracted image exists but the mask is missing: {brain}"

    command = ["bet", str(magnitude), str(prefix), "-m", "-R", "-f", str(fraction)]
    try:
        proc = subprocess.run(command, text=True, capture_output=True, check=False)
    except OSError as exc:
        return "ERROR", str(mask), f"Could not start BET: {exc}"
    if proc.returncode != 0 or not mask.is_file() or mask.stat().st_size == 0:
        detail = (proc.stderr or proc.stdout).strip().replace("\n", " ")
        return "ERROR", str(mask), f"BET failed (rc={proc.returncode}): {detail[:500]}"

    # Verify the generated mask is binary-valued and non-empty using FSL tools.
    try:
        stats = subprocess.run(["fslstats", str(mask), "-R", "-V"],
                               text=True, capture_output=True, check=True)
        values = [float(x) for x in stats.stdout.split()]
        if (len(values) < 4 or values[0] < 0 or values[1] > 1
                or values[1] != 1 or values[2] <= 0):
            return "QC_FAIL", str(mask), f"Unexpected mask range or volume: {stats.stdout.strip()}"
    except (OSError, subprocess.CalledProcessError, ValueError) as exc:
        return "QC_FAIL", str(mask), f"Could not validate the mask: {exc}"
    return "OK", str(mask), f"BET -f {fraction:g}; voxel count={int(values[2])}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worklist", type=Path, default=DEFAULT_WORKLIST)
    parser.add_argument("--sample-root", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--fraction", type=float, default=0.30,
                        help="BET fractional intensity threshold (default: 0.30)")
    parser.add_argument("--overwrite", action="store_true",
                        help="Regenerate and overwrite existing BET outputs")
    parser.add_argument("--dry-run", action="store_true",
                        help="List images to process without running BET")
    args = parser.parse_args()

    if not args.worklist.is_file():
        print(f"ERROR: Worklist does not exist: {args.worklist}", file=sys.stderr)
        return 2
    if not 0 < args.fraction < 1:
        print("ERROR: --fraction must be between 0 and 1", file=sys.stderr)
        return 2
    if shutil.which("bet") is None or shutil.which("fslstats") is None:
        print("ERROR: bet/fslstats not found; load FSL first, for example: module load fsl/6.0.4",
              file=sys.stderr)
        return 2

    try:
        grouped = read_worklist(args.worklist)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    if not args.sample_root.is_dir():
        print(f"ERROR: Sample directory does not exist: {args.sample_root}", file=sys.stderr)
        return 2
    sample_ids = {p.name for p in args.sample_root.iterdir() if p.is_dir()}
    selected_ids = sorted(set(grouped) & sample_ids)
    absent_ids = sorted(set(grouped) - sample_ids)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    manifest = args.sample_root / f"brain_mask_manifest_{stamp}.csv"
    records = []
    counts = defaultdict(int)
    print(f"Worklist image IDs included in the main model and runnable: {len(grouped)}")
    print(f"Image IDs present in sample/ and matched to the worklist: {len(selected_ids)}")
    print(f"Image IDs absent from sample/ and not processed: {len(absent_ids)}")
    print(f"Sample/mask directory: {args.sample_root}")

    for iid in selected_ids:
        magnitude, echo_list, echo_count, echo_basis, error = choose_source(iid, grouped[iid])
        if error:
            status, mask_path, detail = "SKIP", "", error
        elif args.dry_run:
            status = "DRY_RUN"
            mask_path = str(args.sample_root / iid / f"{iid}_M_e1_brain_mask.nii.gz")
            detail = f"First echo={magnitude}; echo count={echo_count}; selection basis={echo_basis}"
        else:
            status, mask_path, detail = run_bet(
                iid, magnitude, args.sample_root, args.fraction, args.overwrite
            )
        counts[status] += 1
        records.append({"image_dir_id": iid, "first_echo_magnitude": str(magnitude or ""),
                        "magnitude_echo_count": echo_count, "magnitude_echo_files": echo_list,
                        "echo_selection_basis": echo_basis, "brain_mask_file": mask_path,
                        "status": status, "detail": detail})
        print(f"{status}\t{iid}\t{detail}")

    with manifest.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("image_dir_id", "first_echo_magnitude",
                                                    "magnitude_echo_count", "magnitude_echo_files",
                                                    "echo_selection_basis", "brain_mask_file",
                                                    "status", "detail"))
        writer.writeheader()
        writer.writerows(records)
    print(f"Summary: {dict(sorted(counts.items()))}")
    print(f"Manifest: {manifest}")
    print("Note: These are draft brain masks generated from magnitude images, not QSM reconstruction-support masks.")
    return 1 if counts["ERROR"] or counts["QC_FAIL"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
