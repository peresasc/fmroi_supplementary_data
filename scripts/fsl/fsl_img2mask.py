#===========================================================================
# FSL Image Threshold ROI Generation (img2mask)
#
#--------------------------------------------------------------------------
# This script generates Region of Interests (ROIs) based on minimum and
# maximum threshold using the FSL software. It assesses FSL's functionality
# and performance in creating ROIs, as reported in the paper "fMROI: a simple
# and adaptable toolbox for easy region-of-interest creation," currently under
# review.
#
# Author: Igor Vaz, 2024.
#===========================================================================

import subprocess
import random
import os
import nibabel as nib
import numpy as np


def get_image_min_max(image_path):
    """Get the minimum and maximum values from a NIfTI image."""
    img = nib.load(image_path)
    data = img.get_fdata()
    return np.min(data), np.max(data)


def generate_random_thresholds(imgmin, imgmax):
    """Generate a pair of random thresholds within the min and max of the image."""
    if imgmin < 0:
        imgmin = 0
    t1 = round(random.uniform(imgmin, imgmax), 4)
    t2 = round(random.uniform(imgmin, imgmax), 4)
    return min(t1, t2), max(t1, t2)


def run_fslmaths(thr, uthr, iteration, outdir, prefix):
    """Run the fslmaths command with the specified thresholds and save the output."""
    # Select the image based on the prefix
    if prefix == "syndata":
        image_path = "./complex-shapes.nii.gz"
    else:
        image_path = "./default_mode_association-test_z_FDR_0.01.nii.gz"
    
    # Adjust the format specifiers for floating-point numbers {uthr:05.1f}
    if prefix == "dmn":
        output_filename = f"fsl-img2mask_srcimg_dmn_threshold_{thr}_{uthr}.nii"
    else:
        output_filename = f"fsl-img2mask_srcimg_{prefix}_threshold_{thr}_{uthr}.nii"
    
    output_path = os.path.join(outdir, output_filename)
    command = f"fslmaths {image_path} -thr {thr} -uthr {uthr} -abs -bin {output_path}"
    subprocess.run(command, shell=True)
    print(f"Iteration {iteration}: fslmaths command run with thr={thr} and uthr={uthr}, output saved to {output_path}")


def main():
    # Run for random thresholds
    outdir_random = "rois_dmn"
    imgmin, imgmax = get_image_min_max("./default_mode_association-test_z_FDR_0.01.nii.gz")
    os.makedirs(outdir_random, exist_ok=True)
    for i in range(100):
        thr, uthr = generate_random_thresholds(imgmin, imgmax)
        run_fslmaths(thr, uthr, i + 1, outdir_random, "dmn")

    # Run for predefined threshold pairs
    outdir_predefined = "rois_syndata"
    os.makedirs(outdir_predefined, exist_ok=True)
    threshold_pairs = [
        (1, 1), (3, 4), (5, 5), (5.2, 5.2), (5.4, 5.4),
        (5.6, 5.6), (5.8, 5.8), (7, 8), (9, 10), (11, 12)
    ]
    for i, (thr, uthr) in enumerate(threshold_pairs, start=1):
        run_fslmaths(thr, uthr, i, outdir_predefined, "syndata")


if __name__ == "__main__":
    main()