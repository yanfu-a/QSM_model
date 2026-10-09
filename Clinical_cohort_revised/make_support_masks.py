#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Write QSM support masks (<IID>_support.nii.gz) and their provenance for P9's --support-dir.

The support mask marks the voxels where the QSM reconstruction produced values. Two kinds
are written, and each mask's JSON record (<IID>_support.json) says which:

* support_type=reconstruction: the reconstruction pipeline's own brain/ROI mask, registered
  with --recon-mask-template. It must already lie on the QSM grid; it is binarized (> 0).
* support_type=estimated_from_qsm: for QSM maps exported without their mask (the cohort's
  vendor QSM). The support is recovered from the QSM image: the nonzero, finite region plus
  zero-valued holes it fully encloses of at most --max-hole-mm3. Small enclosed zeros are
  treated as genuine values (common when QSM is stored as quantized integers); larger or
  open zero regions (masked CSF, outside the brain) stay unsupported. This is an estimate,
  not independent evidence of the reconstruction support.

The record holds the SHA-1 of the source QSM (and reconstruction mask), the method and its
parameters, and the SHA-1 of the written mask. An existing mask is reused only when its
record matches the current QSM, source mask, method and parameters and the mask file is
unchanged; otherwise it is rebuilt. P9 verifies the same record before using a mask.
Masks and records are written to temporary files and renamed, so an interrupted run never
leaves a record describing a partial mask.
"""

import argparse
import datetime
import hashlib
import json
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage

DEFAULT_NATIVE_DIR = Path("/cwStorage/nodecw_group/FY_data/QSM_HUASHAN/native_data")
HERE = Path(__file__).resolve().parent
ESTIMATED_METHOD = "nonzero_plus_enclosed_holes_v1"
RECON_METHOD = "reconstruction_mask_binarized_v1"


def file_sha1(path, chunk=1 << 20):
    """SHA-1 of a file's contents."""
    digest = hashlib.sha1()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def provenance_path(mask_path):
    """<IID>_support.json next to <IID>_support.nii[.gz]."""
    name = Path(mask_path).name
    stem = name[:-7] if name.endswith(".nii.gz") else name[:-4] if name.endswith(".nii") else name
    return Path(mask_path).with_name(stem + ".json")


def write_atomic_bytes(path, data):
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_bytes(data)
    os.replace(tmp, path)


def support_from_qsm(qsm, vox_mm3, max_hole_mm3):
    """Nonzero finite QSM plus enclosed zero holes of at most max_hole_mm3.

    Return (support, n_holes_filled, filled_voxels, largest_filled_hole_mm3).
    """
    nonzero = np.isfinite(qsm) & (qsm != 0)
    holes, n_comp = ndimage.label(~nonzero)       # face-connected zero regions
    if n_comp == 0:
        return nonzero, 0, 0, 0.0
    sizes = np.bincount(holes.ravel(), minlength=n_comp + 1)
    # Regions touching the array border are outside the reconstruction, never holes.
    border = np.unique(np.concatenate([
        holes[0].ravel(), holes[-1].ravel(), holes[:, 0].ravel(), holes[:, -1].ravel(),
        holes[:, :, 0].ravel(), holes[:, :, -1].ravel()]))
    fill = (sizes * vox_mm3 <= max_hole_mm3)
    fill[0] = False
    fill[border] = False
    support = nonzero | fill[holes]
    filled = int(sizes[fill].sum())
    largest = float(sizes[fill].max() * vox_mm3) if fill.any() else 0.0
    return support, int(fill.sum()), filled, largest


def record_matches(record, wanted, mask_path):
    """True when an existing record describes these inputs and the mask file is unchanged."""
    keys = ("support_type", "method", "qsm_sha1", "max_hole_mm3", "source_mask_sha1")
    if any(record.get(k) != wanted.get(k) for k in keys):
        return False
    return Path(mask_path).is_file() and record.get("mask_sha1") == file_sha1(mask_path)


def make_one(oid, qsm_path, out_dir, max_hole_mm3, overwrite, recon_template):
    """Write (or verify and keep) one subject's support mask; return its summary row."""
    out = out_dir / f"{oid}_support.nii.gz"
    prov = provenance_path(out)
    row = {"IID": oid, "qsm_file": qsm_path, "support_mask": str(out)}
    try:
        qsm_sha1 = file_sha1(qsm_path)
        if recon_template:
            source = Path(recon_template.format(IID=oid))
            if not source.is_file():
                return {**row, "status": f"ERROR reconstruction mask missing: {source}"}
            wanted = {"support_type": "reconstruction", "method": RECON_METHOD, "qsm_sha1": qsm_sha1,
                      "max_hole_mm3": None, "source_mask_sha1": file_sha1(source)}
        else:
            source = None
            wanted = {"support_type": "estimated_from_qsm", "method": ESTIMATED_METHOD,
                      "qsm_sha1": qsm_sha1, "max_hole_mm3": float(max_hole_mm3), "source_mask_sha1": None}

        if not overwrite and prov.is_file() and out.is_file():
            record = json.loads(prov.read_text(encoding="utf-8"))
            if record_matches(record, wanted, out):
                return {**row, **{k: record.get(k) for k in SUMMARY_FIELDS}, "status": "kept (provenance verified)"}
            reason = "stale (provenance changed)"
        else:
            reason = "rebuilt" if overwrite and out.is_file() else "written"

        img = nib.load(qsm_path)
        qsm = np.asanyarray(img.dataobj)
        if qsm.ndim > 3:
            qsm = qsm[..., 0]
        nonzero = np.isfinite(qsm) & (qsm != 0)
        vox_mm3 = float(abs(np.linalg.det(img.affine[:3, :3])))
        if source is not None:
            src = nib.load(str(source))
            if src.shape[:3] != img.shape[:3] or not np.allclose(src.affine, img.affine, atol=1e-3):
                return {**row, "status": f"ERROR reconstruction mask is not on the QSM grid: {source}"}
            sdata = np.asanyarray(src.dataobj)
            support = (sdata[..., 0] if sdata.ndim > 3 else sdata) > 0
            n_holes, filled, largest = 0, 0, 0.0
        else:
            support, n_holes, filled, largest = support_from_qsm(qsm, vox_mm3, max_hole_mm3)

        # Same grid as the QSM, so P9's shape/affine check accepts it.
        mask = nib.Nifti1Image(support.astype(np.uint8), img.affine)
        mask.set_sform(img.affine, int(img.header["sform_code"]) or 1)
        mask.set_qform(img.affine, int(img.header["qform_code"]) or 1)
        prov.unlink(missing_ok=True)            # never leave an old record next to a new mask
        tmp = out.with_name(f".{out.name}.tmp.nii.gz")
        nib.save(mask, str(tmp))
        os.replace(tmp, out)
        record = {"IID": oid, **wanted, "qsm_file": str(qsm_path),
                  "source_mask_file": str(source) if source is not None else None,
                  "mask_file": str(out), "mask_sha1": file_sha1(out),
                  "qsm_dtype": str(img.get_data_dtype()),
                  "qsm_integer_valued": bool(np.all(np.mod(qsm[nonzero], 1) == 0)) if nonzero.any() else False,
                  "nonzero_voxels": int(nonzero.sum()), "holes_filled": n_holes, "filled_voxels": filled,
                  "largest_filled_hole_mm3": round(largest, 3), "support_voxels": int(support.sum()),
                  "qsm_nonzero_outside_support": int((nonzero & ~support).sum()),
                  "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
                  "script": "make_support_masks.py"}
        write_atomic_bytes(prov, (json.dumps(record, indent=1) + "\n").encode("utf-8"))
        return {**row, **{k: record.get(k) for k in SUMMARY_FIELDS}, "status": reason}
    except Exception as exc:                                       # noqa: BLE001
        return {**row, "status": f"ERROR {type(exc).__name__}: {exc}"}


SUMMARY_FIELDS = ("support_type", "method", "max_hole_mm3", "qsm_sha1", "mask_sha1", "qsm_dtype",
                  "qsm_integer_valued", "nonzero_voxels", "holes_filled", "filled_voxels",
                  "largest_filled_hole_mm3", "support_voxels", "qsm_nonzero_outside_support")


def main():
    ap = argparse.ArgumentParser(description="Write QSM support masks and provenance for P9 --support-dir")
    ap.add_argument("--worklist", type=Path, default=HERE / "worklist.csv")
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_NATIVE_DIR / "support")
    ap.add_argument("--ids-file", type=Path, default=None,
                    help="Only these subjects (IID or image_dir_id column)")
    ap.add_argument("--max-hole-mm3", type=float, default=10.0,
                    help="Estimated masks: largest enclosed zero region treated as genuine QSM values")
    ap.add_argument("--recon-mask-template", default="",
                    help="Register genuine reconstruction masks instead of estimating, e.g. "
                         "'/data/{IID}/qsm_mask.nii.gz' ({IID} is replaced by image_dir_id)")
    ap.add_argument("--overwrite", action="store_true",
                    help="Rebuild every mask even when its provenance record verifies")
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args()
    if args.recon_mask_template and "{IID}" not in args.recon_mask_template:
        ap.error("--recon-mask-template must contain {IID}")

    wl = pd.read_csv(args.worklist, dtype=str).fillna("")
    if "runnable" in wl.columns:
        wl = wl[wl["runnable"].str.strip().str.lower() == "yes"]
    wl = wl.drop_duplicates("image_dir_id")
    if args.ids_file:
        p = pd.read_csv(args.ids_file, dtype=str).fillna("")
        idcol = next((c for c in ("IID", "image_dir_id") if c in p.columns), None)
        if idcol is None:
            raise SystemExit(f"{args.ids_file} has neither an IID nor an image_dir_id column")
        wl = wl[wl["image_dir_id"].isin(p[idcol])]
    args.out_dir.mkdir(parents=True, exist_ok=True)

    with ThreadPoolExecutor(args.jobs) as ex:
        rows = list(ex.map(lambda t: make_one(t[0], t[1], args.out_dir, args.max_hole_mm3,
                                              args.overwrite, args.recon_mask_template),
                           zip(wl["image_dir_id"], wl["qsm_file"])))
    summary = pd.DataFrame(rows)
    summary_path = args.out_dir / "support_masks_summary.csv"
    if summary_path.is_file():
        # Keep rows of subjects not processed in this run.
        old = pd.read_csv(summary_path, dtype=str)
        summary = pd.concat([old[~old["IID"].isin(summary["IID"])], summary], ignore_index=True)
    write_atomic_bytes(summary_path, summary.to_csv(index=False).encode("utf-8-sig"))

    status = pd.Series([r["status"] for r in rows], dtype=str)
    print(f"Support masks in {args.out_dir}: {int(status.isin(['written', 'rebuilt']).sum())} written, "
          f"{int(status.str.startswith('stale').sum())} rebuilt as stale, "
          f"{int(status.str.startswith('kept').sum())} kept (provenance verified), "
          f"{int(status.str.startswith('ERROR').sum())} failed (summary: {summary_path.name})")
    done = pd.DataFrame([r for r in rows if not r["status"].startswith("ERROR")])
    if len(done) and "filled_voxels" in done:
        frac = done["filled_voxels"].astype(float) / done["nonzero_voxels"].astype(float).clip(lower=1)
        print(f"  enclosed zero voxels added: median fraction {frac.median():.2e}, max {frac.max():.2e}; "
              f"integer-valued QSM: {int(done['qsm_integer_valued'].astype(str).eq('True').sum())}")
    for r in rows:
        if r["status"].startswith("ERROR"):
            print(f"  {r['IID']}: {r['status']}")


if __name__ == "__main__":
    main()
