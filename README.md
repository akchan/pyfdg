# pyfdg
 
SUV calculation from dicom files based on pseudo-code from QIBA SUV subcommittee.

https://qibawiki.rsna.org/index.php/Standardized_Uptake_Value_(SUV)

# Requirements

- numpy
- pydicom

```bash
pip install -r requirements.txt
```

# Usage

## Calculate SUVmax from a dicom file

```
$ python pydicom.py {dicom_file_path}
```

## Use in your python code

First, clone this repository or pyfdg.py file to your project.

```py
import pyfdg

dcm_dir_path = 'path'

# SUV_bw case
vol_suv_bw, voxel_size, matrix_size = pyfdg.read_fdg_pet_dicom_dir(dcm_dir_path)

suv_bw_max = np.max(vol_suv[:])
print('SUV_bw_max:', suv_bw_max)

# SUV_LBM case
patient_height = 1.70  # in meter

vol_suv_lbm, voxel_size, matrix_size = pyfdg.read_fdg_pet_dicom_dir(dcm_dir_path, target='suv_lbm', patient_height=patient_height)

suv_lbm_max = np.max(vol_suv_lbm[:])
print('SUV_bw_max:', suv_lbm_max)

```