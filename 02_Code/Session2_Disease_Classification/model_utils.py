import torch
from torch.utils.data import Dataset
import pandas as pd
import numpy as np
import nibabel as nib
import os
import random

mask = nib.load('/path_to/')
mask_data = mask.get_fdata()
nonzero = mask_data.nonzero()


class Multi_disease_dataset(Dataset):

    def __init__(self, image_dir, image_dir_add, data_csv):
        self.df = pd.read_csv(data_csv,dtype={'sub_id':str})
        self.image_dir = image_dir
        self.image_dir_add = image_dir_add
        self.n_data = len(self.df)

    def __len__(self):
        return self.n_data

    def __getitem__(self, idx):
        row = self.df.iloc[idx]
        image_path = os.path.join(self.image_dir,row['sub_id'],'T1w_brain_MNI_nonlinear.npy')
        if not os.path.exists(image_path): 
            image_path = os.path.join(self.image_dir_add,row['sub_id'],'T1w_brain_MNI_nonlinear.npy')
        if not os.path.exists(image_path):
            image  = nib.load(row['high_res_path']).get_fdata()[nonzero] 
            image_path = os.path.join(self.image_dir_add,row['sub_id'],'T1w_brain_MNI_nonlinear.npy')
            np.save(image_path,image)

        image_1d = np.load(image_path).reshape(-1)
        image_1d = (image_1d - image_1d.mean()) / image_1d.std()
        image = mask_data.copy()
        image[nonzero] = image_1d
        image = image[9:81, 10:100, 1:77]
        image = np.expand_dims(image, axis=0)  # (1, D, H, W)
        label = int(row['label_id'])
        return torch.tensor(image, dtype=torch.float32), torch.tensor(label)



class Multi_disease_dataset_with_clinical(Dataset):

    def __init__(self, image_dir, image_dir_add, recron_image_dir, data_csv, clinical_csv):
        self.df = pd.read_csv(data_csv,dtype={'sub_id':str})
        self.df_cli = pd.read_csv(clinical_csv,dtype={'clinical_sub_id':str})
        self.image_dir = image_dir
        self.recron_image_dir = recron_image_dir
        self.n_data = len(self.df)
        self.n_clinical = len(self.df_cli)
        print(f"Total {self.n_data + self.n_clinical} images for training, containing {self.n_clinical} reconstruction images.")

    def __len__(self):
        return self.n_data + self.n_clinical

    def __getitem__(self, idx):
        if idx < self.n_data:
            row = self.df.iloc[idx]
            image_path = os.path.join(self.image_dir,row['sub_id'],'T1w_brain_MNI_nonlinear.npy')
            if not os.path.exists(image_path): 
                image_path = os.path.join(self.image_dir_add,row['sub_id'],'T1w_brain_MNI_nonlinear.npy')
            if not os.path.exists(image_path):
                image  = nib.load(row['high_res_path']).get_fdata()[nonzero] 
                image_path = os.path.join(self.image_dir_add,row['sub_id'],'T1w_brain_MNI_nonlinear.npy')
                np.save(image_path,image)
        else:
            row = self.df_cli.iloc[idx-self.n_data]
            image_path = os.path.join(self.recron_image_dir,row['clinical_sub_id']+'_'+row['clinical_plane'],'T1W_brain_nonlinear_recon_average_5.npy')

        image_1d = np.load(image_path).reshape(-1)
        image_1d = (image_1d - image_1d.mean()) / image_1d.std()
        image = mask_data.copy()
        image[nonzero] = image_1d
        image = image[9:81, 10:100, 1:77]
        image = np.expand_dims(image, axis=0)  # (1, D, H, W)
        label = int(row['label_id'])
        return torch.tensor(image, dtype=torch.float32), torch.tensor(label)



