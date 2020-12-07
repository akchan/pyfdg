# qiba_fdg_suv

SUV calculation from dicom files based on pseudo-code from QIBA SUV subcommittee.

https://qibawiki.rsna.org/index.php/Standardized_Uptake_Value_(SUV)

# Requirements

- numpy
- pydicom

# Usage

## Calculate SUVmax

```
$ python suv.py {dicom_file_path}
```

## Use in your python code

First, clone this repository or suv.py file to your project.

```py
import suv

import pydicom

dcm_path = 'path_to_dicom'

dcm = pydicom.read_file(dcm_path)

img_suv = suv.calc_suv_qiba(dcm.pixel_array, dcm)

suv_max = np.max(img_suv[:])
print('SUVmax:', suv_max)

```