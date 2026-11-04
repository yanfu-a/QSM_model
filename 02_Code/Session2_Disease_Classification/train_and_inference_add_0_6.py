import os
import random
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from sklearn.metrics import accuracy_score, roc_auc_score, balanced_accuracy_score
from resnet import resnet10 as resnet_model
from model_utils import Multi_disease_dataset, Multi_disease_dataset_with_clinical

def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def create_file_add_ratio(save_dir, file_path, add_ratio, seed):
    data_file = pd.read_csv(file_path, dtype={'clinical_sub_id': str})
    select_num = int(add_ratio * len(data_file))
    sub_data_file = data_file.sample(n=select_num, random_state=seed)

    new_save_dir = os.path.join(save_dir, 'add_ratio_' + str(add_ratio))
    os.makedirs(new_save_dir, exist_ok=True)

    sub_file_path = os.path.join(new_save_dir, 'sub_file.csv')
    sub_data_file.to_csv(sub_file_path, index=False)
    return new_save_dir, sub_file_path


def train_epoch(model, loader, optimizer, criterion, device):
    model.train()
    running_loss = 0.0

    for images, targets in loader:
        images = images.to(device)
        targets = targets.to(device).long()

        optimizer.zero_grad()
        outputs = model(images)
        loss = criterion(outputs, targets)
        loss.backward()
        optimizer.step()
        running_loss += loss.item() * images.size(0)

    return running_loss / len(loader.dataset)


def eval_epoch(model, loader, criterion, device, n_classes, return_predictions=False):
    model.eval()
    all_preds, all_targets, all_probs = [], [], []
    running_loss = 0.0

    with torch.no_grad():
        for images, labels in loader:
            images = images.to(device)
            labels = labels.to(device).long()

            outputs = model(images)
            loss = criterion(outputs, labels)
            probs = torch.softmax(outputs, dim=1)
            preds = torch.argmax(probs, dim=1)

            running_loss += loss.item() * images.size(0)
            all_preds.extend(preds.cpu().numpy())
            all_targets.extend(labels.cpu().numpy())
            all_probs.extend(probs.cpu().numpy())

    all_targets = np.array(all_targets)
    all_preds = np.array(all_preds)
    all_probs = np.array(all_probs)

    if len(all_targets) == 0:
        raise ValueError("Evaluation dataset is empty.")

    if all_probs.shape[1] != n_classes or np.any((all_targets < 0) | (all_targets >= n_classes)):
        raise ValueError("Model outputs or labels do not match n_classes.")

    acc = accuracy_score(all_targets, all_preds)
    balanced_acc = balanced_accuracy_score(all_targets, all_preds)
    val_loss = running_loss / len(loader.dataset)

    # Macro OvR AUROC: only evaluate classes present in the dataset
    present_classes = np.unique(all_targets)

    if len(present_classes) >= 2:
        auc_list = []
        for c in present_classes:
            binary_targets = (all_targets == c).astype(int)
            auc_c = roc_auc_score(binary_targets, all_probs[:, c])
            auc_list.append(auc_c)
        auc = float(np.mean(auc_list))
    else:
        auc = np.nan

    if return_predictions:
        return acc, auc, balanced_acc, val_loss, all_preds, all_probs
    else:
        return acc, auc, balanced_acc, val_loss


def main(image_dir, image_dir_add, recron_image_dir, train_csv, clinical_csv, val_csv, test_csv,
         model_save_dir, n_classes=12, batch_size=16, lr=1e-4, epochs=50, device='cuda', seed=42):

    set_seed(seed)
    os.makedirs(model_save_dir, exist_ok=True)

    test_data = pd.read_csv(test_csv)
    results = []
    results_path = os.path.join(model_save_dir, 'results.csv')

    for ratio_idx in range(2, 11):
        add_ratio = ratio_idx / 10
        set_seed(seed)

        save_dir, sub_clinical_csv = create_file_add_ratio(model_save_dir, clinical_csv, add_ratio, seed)

        weight_path = os.path.join(save_dir, 'best_model.pth')
        training_log_path = os.path.join(save_dir, 'training_log.csv')
        test_results_path = os.path.join(save_dir, 'test_results.csv')

        train_ds = Multi_disease_dataset_with_clinical(image_dir, image_dir_add, recron_image_dir, train_csv, sub_clinical_csv)
        val_ds = Multi_disease_dataset(image_dir, image_dir_add, val_csv)
        test_ds = Multi_disease_dataset(image_dir, image_dir_add, test_csv)

        train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True, num_workers=8)
        val_loader = DataLoader(val_ds, batch_size=batch_size, shuffle=False, num_workers=8)
        test_loader = DataLoader(test_ds, batch_size=batch_size, shuffle=False, num_workers=8)

        model = resnet_model(inchan=1, num_classes=n_classes).to(device)
        criterion = nn.CrossEntropyLoss()
        optimizer = optim.Adam(model.parameters(), lr=lr)

        # Use AUROC when >= 2 classes are present; otherwise use validation loss
        best_score = -np.inf
        best_val_auc = np.nan
        best_val_loss = np.nan
        best_epoch = 0
        training_log = []

        for epoch in range(1, epochs + 1):
            train_loss = train_epoch(model, train_loader, optimizer, criterion, device)
            val_acc, val_auc, val_b_acc, val_loss = eval_epoch(model, val_loader, criterion, device, n_classes)

            training_log.append({
                'epoch': epoch,
                'train_loss': train_loss,
                'val_loss': val_loss,
                'val_acc': val_acc,
                'val_auc': val_auc,
                'val_b_acc': val_b_acc
            })

            print(f"Epoch {epoch}/{epochs}: Train Loss={train_loss:.4f} | Val Loss={val_loss:.4f} | "
                  f"Val Acc={val_acc:.4f} | Val AUC={val_auc:.4f} | Val B_Acc={val_b_acc:.4f}")

            # Select best model
            score = val_auc if np.isfinite(val_auc) else -val_loss

            if not np.isfinite(score):
                raise ValueError(f"Invalid validation metric at epoch {epoch}")

            if score > best_score:
                best_score = score
                best_val_auc = val_auc
                best_val_loss = val_loss
                best_epoch = epoch
                torch.save(model.state_dict(), weight_path)
                print(f">>> New best model saved at {weight_path}")

        pd.DataFrame(training_log).to_csv(training_log_path, index=False)

        # Test
        model.load_state_dict(torch.load(weight_path, map_location=device, weights_only=True))
        test_acc, test_auc, test_b_acc, test_loss, test_preds, test_probs = eval_epoch(
            model, test_loader, criterion, device, n_classes, return_predictions=True
        )

        if len(test_data) != len(test_preds):
            raise ValueError("Test CSV and predictions have different lengths.")

        test_output = test_data.copy()
        test_output['pred_label'] = test_preds

        for c in range(n_classes):
            test_output[f'prob_class_{c}'] = test_probs[:, c]

        test_output.to_csv(test_results_path, index=False)

        results.append({
            'add_ratio': add_ratio,
            'best_epoch': best_epoch,
            'best_val_auc': best_val_auc,
            'best_val_loss': best_val_loss,
            'test_acc': test_acc,
            'test_auc': test_auc,
            'test_b_acc': test_b_acc,
            'test_loss': test_loss
        })

        pd.DataFrame(results).to_csv(results_path, index=False)
        print(f"Ratio={add_ratio} | Test ACC={test_acc:.4f} | Test AUC={test_auc:.4f} | Test B_Acc={test_b_acc:.4f}")


if __name__ == '__main__':
    config = {
        'image_dir': '/path_to/',
        'image_dir_add': '/path_to/',
        'recron_image_dir': '/path_to/',
        'train_csv': '/path_to/',
        'clinical_csv': '/path_to/',
        'val_csv': '/path_to/',
        'test_csv': '/path_to/',
        'model_save_dir': '/path_to/',
        'n_classes': 7,
        'batch_size': 16,
        'lr': 1e-4,
        'epochs': 50,
        'device': 'cuda' if torch.cuda.is_available() else 'cpu',
        'seed': 42
    }

    main(**config)