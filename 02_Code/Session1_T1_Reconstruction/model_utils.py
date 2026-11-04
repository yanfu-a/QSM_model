# This file contains modulues used for Dit_1d
import os
import torch
import torch.nn as nn
import numpy as np
import pandas as pd
from torch.utils.data import Dataset
import nibabel as nib
import torch.nn.functional as F

mask = nib.load('/cpfs01/projects-HDD/cfff-afe2df89e32e_HDD/public/work_files_yc/BrainGen/yc_code/preprocess/MNI152_T1_2mm_brain.nii.gz')
mask_data = mask.get_fdata()
nonzero = mask_data.nonzero()

class PatchEmbed_1d(nn.Module):
    """ 1D Image to Patch Embedding
    """
    def __init__(self, input_size, patch_size=256, padding_size=155, in_chans=1, embed_dim=768, norm_layer=None):
        super().__init__()

        self.input_size = input_size
        self.patch_size = patch_size
        self.num_patches = (self.input_size + padding_size) // self.patch_size
        self.padding_size = padding_size

        assert self.num_patches * self.patch_size == self.input_size + padding_size, \
            f"Hidden size error: input size extend ({self.input_size + padding_size}) must be divisible by hidden size ({self.patch_size})"

        self.proj = nn.Conv1d(in_chans, embed_dim, kernel_size=patch_size, stride=patch_size)
        self.norm = norm_layer(embed_dim) if norm_layer else nn.Identity()

    def forward(self, x):
        x = F.pad(x, (0, self.padding_size))
        x = self.proj(x)
        x = x.transpose(1, 2)  
        x = self.norm(x)
        return x


class UKBDataset_T1_sr(Dataset):
    def __init__(self, data_path, data_add_path, data_file,transform=None):

        # Read the CSV file to get labels
        self.data_path = data_path
        self.data_add_path = data_add_path
        self.data_file = pd.read_csv(data_file,dtype={'clinical_sub_id': str})
        self.transform = transform

    def __len__(self):
        return len(self.data_file)
    
    def __getitem__(self, idx):

        # Get the row for the current index
        row = self.data_file.iloc[idx]
        eid = row['clinical_sub_id']
        modality = 11    # target modality
        
        age = 0
        sex = 2  

        # T1 image: 1mm T1 image
        # T1 image clinical : clinical T1w image

        
        T1_image_path = os.path.join(self.data_path, str(eid), f"T1w_brain_MNI_nonlinear.npy")
        if not os.path.exists(T1_image_path):
            T1_image_path = os.path.join(self.data_add_path, str(eid), f"T1w_brain_MNI_nonlinear.npy")
        if not os.path.exists(T1_image_path):
            os.makedirs(os.path.join(self.data_add_path, str(eid)))
            T1_image_path = os.path.join(self.data_add_path, str(eid), f"T1w_brain_MNI_nonlinear.npy")
            T1_image = nib.load(row['high_res_path']).get_fdata()[nonzero]
            np.save(T1_image_path, T1_image)

        T1_image = np.load(T1_image_path)
        T1_image = np.nan_to_num(T1_image, nan=0.0)
        T1_image = torch.tensor(T1_image, dtype=torch.float32)


        T1_image_clinical_path = os.path.join(self.data_path, str(eid), row['clinical_plane']+f"_T1W_brain_nonlinear.npy")
        if not os.path.exists(T1_image_clinical_path):
            T1_image_clinical_path = os.path.join(self.data_add_path, str(eid), row['clinical_plane']+f"_T1W_brain_nonlinear.npy")
        if not os.path.exists(T1_image_clinical_path):
            os.makedirs(os.path.join(self.data_add_path, str(eid)))
            T1_image_clinical_path = os.path.join(self.data_add_path, str(eid), row['clinical_plane']+f"_T1W_brain_nonlinear.npy")
            T1_image_clinical = nib.load(row['clinical_path']).get_fdata()[nonzero]
            np.save(T1_image_clinical_path, T1_image_clinical)

        T1_image_clinical = np.load(T1_image_clinical_path)
        T1_image_clinical = np.nan_to_num(T1_image_clinical, nan=0.0)
        T1_image_clinical = torch.tensor(T1_image_clinical, dtype=torch.float32)

        if self.transform:
            T1_image = (T1_image - T1_image.mean())/T1_image.std()
            T1_image_clinical = (T1_image_clinical - T1_image_clinical.mean())/T1_image_clinical.std()

        # Convert to tensor
        data = {}
        age = torch.tensor(age, dtype=torch.float32)
        sex = torch.tensor(sex, dtype=torch.int32)
        modality = torch.tensor(modality,dtype=torch.int32)

        data ={ 'T1_image':T1_image.unsqueeze(0), 'T1_image_clinical':T1_image_clinical.unsqueeze(0), 
                'age': age,'sex': sex,'mod':modality
               }
        return data 


class UKBDataset_T1_sr_eval(Dataset):
    def __init__(self, data_path, data_add_path, data_file,transform=None):

        # Read the CSV file to get labels
        self.data_path = data_path
        self.data_add_path = data_add_path
        self.data_file = pd.read_csv(data_file,dtype={'clinical_sub_id': str})
        self.transform = transform

    def __len__(self):
        return len(self.data_file)
    
    def __getitem__(self, idx):

        # Get the row for the current index
        row = self.data_file.iloc[idx]
        eid = row['clinical_sub_id']

        modality = 11    # target modality
        
        age = 0
        sex = 2  

        # T1 image: 1mm T1 image
        # T1 image clinical : clinical T1w image

        T1_image_clinical_path = os.path.join(self.data_path, str(eid), row['clinical_plane']+f"_T1W_brain_nonlinear.npy")

        T1_image_clinical_path = os.path.join(self.data_path, str(eid), row['clinical_plane']+f"_T1W_brain_nonlinear.npy")
        if not os.path.exists(T1_image_clinical_path):
            T1_image_clinical_path = os.path.join(self.data_add_path, str(eid), row['clinical_plane']+f"_T1W_brain_nonlinear.npy")
        if not os.path.exists(T1_image_clinical_path):
            os.makedirs(os.path.join(self.data_add_path, str(eid)))
            T1_image_clinical_path = os.path.join(self.data_add_path, str(eid), row['clinical_plane']+f"_T1W_brain_nonlinear.npy")
            T1_image_clinical = nib.load(row['clinical_path']).get_fdata()[nonzero]
            np.save(T1_image_clinical_path, T1_image_clinical)

        T1_image_clinical = np.load(T1_image_clinical_path)
        T1_image_clinical = np.nan_to_num(T1_image_clinical, nan=0.0)

        if self.transform:
            T1_image_clinical = (T1_image_clinical - T1_image_clinical.mean())/T1_image_clinical.std()
            T1_image_clinical = torch.tensor(T1_image_clinical, dtype=torch.float32).unsqueeze(0)
        
        # Convert to tensor
        data = {}
        age = torch.tensor(age, dtype=torch.float32)
        sex = torch.tensor(sex, dtype=torch.int32)
        modality = torch.tensor(modality,dtype=torch.int32)

        data ={ 'eid':eid,'plane':row['clinical_plane'], 'T1_image_clinical':T1_image_clinical, 
                'age': age,'sex': sex,'mod':modality,
               }

        return data 