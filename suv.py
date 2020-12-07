#!/usr/bin/env python
# coding: UTF-8

import sys

import numpy as np
import pydicom


def dcm_tm_to_sec(tm_str):
    sec = int(tm_str[0:2]) * 60 * 60 \
        + int(tm_str[2:4]) * 60 \
        + int(tm_str[4:6])

    if "." == tm_str[6]:
        sec += float(tm_str[6:13])

    return sec


def dcm_dt_to_time_sec(dt_str):
    return dcm_tm_to_sec(dt_str[8:])


def calc_suv_qiba(img, dcm):
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

                # Radiopharmaceutical start datetime [0x0018, 0x1078]
                # radiopharmaceutical information sequence [0x0054,0x0016]
                if [0x0018, 0x1078] in dcm[0x0054, 0x0016][0]:
                    start_time = dcm_dt_to_time_sec(dcm[0x0054, 0x0016][0][0x0018, 0x1078].value)
                else:
                    # Date is not explicit in this case.
                    start_time = dcm_dt_to_time_sec(dcm[0x0054, 0x0016][0][0x0018, 0x1072].value)

                decay_time = scan_time - start_time

                # Radionuclide Total Dose is NOT corrected for residual dose in syringe, which is ignored here.
                injected_dose = float(dcm[0x0054, 0x0016][0][0x0018, 0x1074].value)

                decayed_dose = injected_dose * np.power(2, -1 * decay_time / half_life)

                weight = float(dcm[0x0010, 0x1030].value)

                suv_bw_scale_factor = weight * 1000 / decayed_dose

                vol_rescaled = pydicom.pixel_data_handlers.apply_modality_lut(vol, dcm)

                suv_bw = vol_rescaled * suv_bw_scale_factor

                # Rescale Intercept is required to be 0 for PET, but use it just in case
                # Rescale slope may vary per slice (GE), and cannot be assumed to be constant for the entire volume

            # This may be post-processed
            else:
                raise NotImplementedError('scan time detection for a post-processed dicom file is not implemented!')

        else:
            raise NotImplementedError('Not inmplemented except for unit of BQML')
    else:
        raise AttributeError('This dicom seems not to be of FDG-PET.')

    return suv_bw


def main(dcm_path):
    dcm = pydicom.read_file(dcm_path)

    img_raw = dcm.pixel_array

    img_suv = calc_suv_qiba(img_raw, dcm)

    print('SUVmax:', np.max(img_suv[:]))


if __name__ == '__main__':
    dcm_path = sys.argv[1]
    main(dcm_path)
