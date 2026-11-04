"""Seven-class ResNet10 inference for UKB T1_brain_to_MNI images.

Example: python inference_ukb_0_6.py --batch-size 16
The output CSV is append-only and an interrupted run can be resumed.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import torch
from scipy.ndimage import affine_transform
from tqdm import tqdm

from resnet import resnet10


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_EIDS = PROJECT_ROOT / "labels/ukb_labels/ukb_t1_sub.csv"
DEFAULT_LABELS = PROJECT_ROOT / "labels/label_id_0_6/classification_train.csv"
DEFAULT_IMAGES = Path(
    "/cpfs01/projects-HDD/cfff-afe2df89e32e_HDD/public/gongwk/ukb_multimodal_new"
)
DEFAULT_MASK = Path(
    "/cpfs01/projects-HDD/cfff-afe2df89e32e_HDD/public/work_files_yc/"
    "BrainGen/yc_code/preprocess/MNI152_T1_2mm_brain.nii.gz"
)
DEFAULT_WEIGHTS = (
    Path(__file__).resolve().parent
    / "results_add_step_40k_lr_1e-4_label_id_0_6/add_ratio_1.0/best_model.pth"
)
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "ukb_inference/ukb_probs_add_ratio_1.0.csv"
N_CLASSES = 7
FAILURE_COLUMNS = ["eid", "image_path", "status", "error"]
CROP = (slice(9, 81), slice(10, 100), slice(1, 77))


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eids", type=Path, default=DEFAULT_EIDS)
    parser.add_argument("--labels", type=Path, default=DEFAULT_LABELS)
    parser.add_argument("--image-dir", type=Path, default=DEFAULT_IMAGES)
    parser.add_argument("--mask", type=Path, default=DEFAULT_MASK)
    parser.add_argument("--weights", type=Path, default=DEFAULT_WEIGHTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--limit", type=int, help="Process only the first N eids for a smoke test")
    return parser.parse_args()


def checkpoint_hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_mask(path):
    image = nib.load(str(path))
    mask = np.asarray(image.dataobj) != 0
    if mask.shape != (91, 109, 91) or not mask.any():
        raise ValueError(f"Unexpected training mask shape or empty mask: {mask.shape}")
    return mask, image.affine


def preprocess(image_path, mask, mask_affine):
    """Match model_utils.Multi_disease_dataset: masked z-score then crop."""
    image = nib.load(str(image_path))
    if len(image.shape) != 3:
        raise ValueError(f"Expected a 3D image, got shape {image.shape}")

    data = np.asarray(image.dataobj, dtype=np.float32)
    target_shape = mask.shape
    if data.shape == target_shape and np.allclose(image.affine, mask_affine, atol=1e-4):
        aligned = data
    elif (data.shape == tuple(2 * n for n in target_shape)
          and np.allclose(image.affine @ np.diag([2, 2, 2, 1]), mask_affine, atol=1e-4)):
        # The usual UKB 1 mm grid lands exactly on the training 2 mm grid.
        aligned = data[::2, ::2, ::2]
    else:
        transform = np.linalg.inv(image.affine) @ mask_affine
        aligned = affine_transform(
            data, transform[:3, :3], offset=transform[:3, 3],
            output_shape=target_shape, order=1, mode="constant", cval=0,
            prefilter=False, output=np.float32,
        )

    values = aligned[mask]
    mean = values.mean(dtype=np.float64)
    std = values.std(dtype=np.float64)
    if not np.isfinite(std) or std <= 0:
        raise ValueError("Masked image intensities have zero or invalid variance")

    volume = np.zeros(target_shape, dtype=np.float32)
    volume[mask] = (values - mean) / std
    cropped = volume[CROP]
    if cropped.shape != (72, 90, 76) or not np.isfinite(cropped).all():
        raise ValueError("Preprocessed image has an invalid shape or non-finite values")
    return torch.from_numpy(np.ascontiguousarray(cropped[None]))


def load_eids(path, limit):
    frame = pd.read_csv(path, dtype={"eid": str})
    if "eid" not in frame.columns:
        raise ValueError(f"Missing eid column in {path}")
    eids = frame["eid"].str.strip()
    if eids.isna().any() or eids.eq("").any() or eids.duplicated().any():
        raise ValueError("Eid column contains an empty or duplicate value")
    if limit is not None:
        if limit < 1:
            raise ValueError("--limit must be positive")
        eids = eids.iloc[:limit]
    return eids.tolist()


def load_class_names(path):
    frame = pd.read_csv(path, usecols=["label_id", "disease_label"])
    pairs = frame.drop_duplicates().sort_values("label_id")
    if (pairs["label_id"].tolist() != list(range(N_CLASSES))
            or pairs["disease_label"].isna().any()
            or pairs["disease_label"].duplicated().any()):
        raise ValueError(f"Expected one unique disease_label for each label_id 0–6 in {path}")
    return pairs["disease_label"].astype(str).tolist()


def completed_eids(path, columns):
    if not path.exists() or path.stat().st_size == 0:
        return set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != columns:
            raise ValueError(f"Existing output has a different schema: {path}")
        return {row["eid"] for row in reader}


def predict_chunk(eids, image_dir, mask, mask_affine, model, device, class_names):
    valid_eids = []
    tensors = []
    failures = []
    for eid in eids:
        image_path = image_dir / eid / "T1_brain_to_MNI.nii.gz"
        if not image_path.is_file():
            failures.append(dict(eid=eid, image_path=str(image_path), status="missing", error="Image file not found"))
        else:
            try:
                tensor = preprocess(image_path, mask, mask_affine)
                valid_eids.append(eid)
                tensors.append(tensor)
            except (OSError, ValueError, RuntimeError) as exc:
                failures.append(dict(eid=eid, image_path=str(image_path), status="error", error=f"{type(exc).__name__}: {exc}"))

    predictions = []
    if tensors:
        with torch.inference_mode():
            logits = model(torch.stack(tensors).to(device))
            if logits.shape != (len(tensors), N_CLASSES):
                raise ValueError(f"Unexpected model output shape: {tuple(logits.shape)}")
            probs = torch.softmax(logits, dim=1).cpu().numpy()
        if not np.isfinite(probs).all():
            raise ValueError("Model produced non-finite probabilities")
        for eid, prob in zip(valid_eids, probs):
            predictions.append({"eid": eid, **{
                name: float(value) for name, value in zip(class_names, prob)
            }})
    return predictions, failures


def main():
    args = parse_args()
    if args.batch_size < 1:
        raise ValueError("--batch-size must be positive")
    eids = load_eids(args.eids, args.limit)
    class_names = load_class_names(args.labels)
    output_columns = ["eid", *class_names]
    weight_hash = checkpoint_hash(args.weights)
    meta_path = args.output.with_suffix(".meta.json")
    metadata = {"checkpoint_sha256": weight_hash, "class_names": class_names}
    if meta_path.exists():
        if json.loads(meta_path.read_text()) != metadata:
            raise ValueError(f"Existing output metadata differs from this run: {meta_path}")
    elif args.output.exists() and args.output.stat().st_size:
        raise ValueError(f"Existing output has no checkpoint metadata: {args.output}")
    completed = completed_eids(args.output, output_columns)
    pending = [eid for eid in eids if eid not in completed]
    print(f"EIDs: {len(eids)}, already written: {len(eids) - len(pending)}, pending: {len(pending)}", flush=True)
    if not pending:
        return

    mask, mask_affine = load_mask(args.mask)
    device = torch.device(args.device)
    model = resnet10(inchan=1, num_classes=N_CLASSES)
    state_dict = torch.load(args.weights, map_location="cpu", weights_only=True)
    if state_dict["fc.weight"].shape[0] != N_CLASSES:
        raise ValueError(f"Checkpoint has {state_dict['fc.weight'].shape[0]} classes; expected {N_CLASSES}: {args.weights}")
    model.load_state_dict(state_dict)
    model.to(device).eval()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    meta_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n")
    failures_path = args.output.with_name(args.output.stem + "_failures_v2.csv")
    write_header = not args.output.exists() or args.output.stat().st_size == 0
    write_failures_header = not failures_path.exists() or failures_path.stat().st_size == 0
    logged_failures = completed_eids(failures_path, FAILURE_COLUMNS)
    with args.output.open("a", newline="") as handle, failures_path.open("a", newline="") as failures_handle:
        writer = csv.DictWriter(handle, fieldnames=output_columns)
        failures_writer = csv.DictWriter(failures_handle, fieldnames=FAILURE_COLUMNS)
        if write_header:
            writer.writeheader()
        if write_failures_header:
            failures_writer.writeheader()
        with tqdm(total=len(pending), desc="UKB inference", unit="image") as progress:
            for start in range(0, len(pending), args.batch_size):
                chunk = pending[start:start + args.batch_size]
                predictions, failures = predict_chunk(chunk, args.image_dir, mask, mask_affine, model, device, class_names)
                writer.writerows(predictions)
                failures_writer.writerows(row for row in failures if row["eid"] not in logged_failures)
                logged_failures.update(row["eid"] for row in failures)
                handle.flush()
                failures_handle.flush()
                progress.update(len(chunk))
    print(f"Saved: {args.output}")
    print(f"Missing/failed images: {failures_path}")


if __name__ == "__main__":
    main()