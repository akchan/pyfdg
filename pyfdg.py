#!/usr/bin/env python
# coding: UTF-8


import glob
import os
import sys

import numpy as np
import pydicom


def list_dicom_files(dir_path, recursive=False):
    # return [str(p) for p in Path(dir_path).glob('*', recursive=True) if os.path.isfile(p) and pydicom.misc.is_dicom(str(p))]
    dir_path = os.path.join(dir_path, '**/*')
    return [str(p) for p in glob.glob(dir_path, recursive=recursive) if os.path.isfile(p) and pydicom.misc.is_dicom(str(p))]


def sort_dcm_files_along_slice_loc(dcm_path_list, omit_slice_zip=True):
    dcm_path_list = [str(p) for p in dcm_path_list]

    assert np.all([pydicom.misc.is_dicom(p) for p in dcm_path_list]), 'All given path should be dicom file.'

    def read_slice_loc(dcm_path):
        dcm = pydicom.read_file(dcm_path)
        loc_z = float(dcm.SliceLocation)
        # loc_z = float(dcm.ImagePositionPatient[2])

        return loc_z

    sorted_dcm_path_list = sorted(dcm_path_list, key=read_slice_loc)

    if omit_slice_zip:
        # Read a file to get scan parameters
        dcm = pydicom.read_file(str(dcm_path_list[0]))

        voxel_size = np.array([float(s) for s in (*dcm.PixelSpacing, dcm.SliceThickness)])

        slice_loc_min = pydicom.read_file(dcm_path_list[0]).SliceLocation
        slice_loc_max = pydicom.read_file(dcm_path_list[-1]).SliceLocation

        # Ignore ZIP (zero fill interpolation along z axis)
        n_slices_ideal = int(
            np.round((slice_loc_max - slice_loc_min) / voxel_size[2])) + 1

        nz = len(dcm_path_list)

        zip_factor = nz / n_slices_ideal

        dcm_path_list = dcm_path_list[::int(zip_factor)]

    return sorted_dcm_path_list


def read_fdg_pet_dicom(dcm_path, target='suv_bw', patient_height=None):
    '''
    Parameters
    ==========

    - target: one of these, suv, sul, raw
    '''
    dcm = pydicom.read_file(dcm_path)

    if target == 'suv_bw' or target == 'suv_lbm':
        # [0x0010, 0x1030] -> patient weight in kg
        body_weight = float(dcm[0x0010, 0x1030].value)

        administered_dose = calc_administered_dose_with_decay(dcm)

        if target == 'suv_bw':
            weight = body_weight
        elif target == 'suv_lbm':
            assert patient_height is not None, 'patient_height is required for SUL.'
            sex = dcm.PatientSex

            if sex == 'M':
                sex = 1
            elif sex == 'F':
                sex = 0
            else:
                raise AttributeError('Invalid sex type:{}'.format(sex))

            weight = calc_lbm(patient_height, body_weight, sex)

        # in g
        weight_g = weight * 1000

        # g/Bq
        suv_scale_factor = weight_g / administered_dose

        # in Bq/mL (Confirmation can be made by checking [0054,1001] Units)
        img_raw = dcm.pixel_array

        # Rescale Intercept is required to be 0 for PET, but use it just in case
        # Rescale slope may vary per slice (GE), and cannot be assumed to be constant for the entire volume
        img_rescaled = pydicom.pixel_data_handlers.apply_modality_lut(img_raw, dcm)

        suv = img_rescaled * suv_scale_factor

        img = suv

    elif target == 'raw':
        img = dcm.pixel_array
    else:
        raise NotImplementedError('Invalid target type:{}'.format(target))

    return img


def calc_lbm(patient_height_m, body_weight_kg, sex):
    '''
    Parameters
    ==========

    - patient_height_m: in m
    - body_weight_kg: in kg
    - sex: 0 for female, 1 for male

    Return
    ======

    lean body mass in kg
    '''
    bmi = body_weight_kg / (patient_height_m * patient_height_m)

    if sex == 0:  # female
        lbm = 9270 * body_weight_kg / (6680 + 216 * bmi)
    elif sex == 1:  # male
        lbm = 9270 * body_weight_kg / (8780 + 244 * bmi)

    return lbm


def read_fdg_pet_dicom_dir(dir_path, recursive=False,
                           target='suv_bw',
                           patient_height=None,
                           omit_slice_zip2=True,
                           dtype=np.float,
                           verbose=False):
    '''
    Parameters
    ==========

    - target: one of these, suv_bw, suv_lbm, raw
    - patient_height: in m
    '''
    # Ignore dot file (i.e. invisible file)
    # paths = list(Path(dir_path).glob('[!.]*.dcm'))
    dcm_path_list = list_dicom_files(dir_path, recursive=recursive)

    if 0 == len(dcm_path_list):
        raise ValueError('An empty directory was passed. {}'.format(dir_path))

    dcm_path_list = sort_dcm_files_along_slice_loc(dcm_path_list,
                                                   omit_slice_zip2=omit_slice_zip2)

    # Read a file to get scan parameters
    dcm = pydicom.read_file(str(dcm_path_list[0]))
    matrix_size = np.array([dcm.Columns, dcm.Rows, len(dcm_path_list)])
    voxel_size = np.array([float(s)
                           for s in (*dcm.PixelSpacing, dcm.SliceThickness)])

    if verbose:
        print('matrix_size:', matrix_size)
        print('voxel_size:', voxel_size)

    vol = np.zeros(matrix_size, dtype=dtype)

    for i, dcm_path in enumerate(dcm_path_list):
        vol[:,:,i] = read_fdg_pet_dicom(dcm_path, target=target, patient_height=patient_height)

    return vol, voxel_size, matrix_size


def dcm_tm_to_sec(tm_str):
    sec = int(tm_str[0:2]) * 60 * 60 \
        + int(tm_str[2:4]) * 60 \
        + int(tm_str[4:6])

    if len(tm_str) >= 7 and "." == tm_str[6]:
        sec += float(tm_str[6:13])

    return sec


def dcm_dt_to_time_sec(dt_str):
    return dcm_tm_to_sec(dt_str[8:])


def calc_administered_dose_with_decay(dcm):
    '''
    Description
    ===========

    Based on QIBA SUV technical subcommiteee

    https://qibawiki.rsna.org/index.php/Standardized_Uptake_Value_(SUV)

    This is an implementation for 'happy path only'.
    '''
    series_date = dcm[0x0008, 0x0021].value
    series_time = dcm[0x0008, 0x0031].value
    acquisition_date = dcm[0x0008, 0x0022].value
    acquisition_time = dcm[0x0008, 0x0032].value

    corrected_image = dcm[0x0028, 0x0051].value
    decay_correction = dcm[0x0054, 0x1102].value
    units = dcm[0x0054, 0x1001].value

    # Radionuclide half life in radiopharmaceutical information sequence (in seconds)
    half_life = float(dcm[0x0054, 0x0016][0][0x0018, 0x1075].value)

    if 'ATTN' in corrected_image and \
       'DECY' in corrected_image and \
       'START' == decay_correction:

        if 'BQML' == units:

            # Be careful. These are comparison of strings.
            if series_date <= acquisition_date and \
               series_time <= acquisition_time:
                scan_time = dcm_tm_to_sec(series_time)

                # [0x0018, 0x1078] -> Radiopharmaceutical start datetime
                # [0x0054, 0x0016] -> Radiopharmaceutical information sequence
                if [0x0018, 0x1078] in dcm[0x0054, 0x0016][0]:
                    # in sec
                    start_time = dcm_dt_to_time_sec(dcm[0x0054, 0x0016][0][0x0018, 0x1078].value)
                else:
                    # Date is not explicit in this case.
                    # in sec
                    start_time = dcm_dt_to_time_sec(dcm[0x0054, 0x0016][0][0x0018, 0x1072].value)

                # in sec
                decay_time = scan_time - start_time

                # Radionuclide Total Dose is NOT corrected for residual dose in syringe, which is ignored here.
                # in Bq
                injected_dose = float(dcm[0x0054, 0x0016][0][0x0018, 0x1074].value)

                # in Bq
                decayed_dose = injected_dose * np.power(2, -1 * decay_time / half_life)

            # This may be post-processed
            else:
                raise NotImplementedError('scan time detection for a post-processed dicom file is not implemented!')

        else:
            raise NotImplementedError('Not inmplemented except for unit of BQML')
    else:
        raise AttributeError('This dicom seems not to be of FDG-PET.')

    return decayed_dose


if __name__ == '__main__':
    dcm_path = sys.argv[1]

    img_suv = read_fdg_pet_dicom(dcm_path)

    suv_max = np.max(img_suv)
    print('suv_max', suv_max)
