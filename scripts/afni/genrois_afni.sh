#bin/bash!

#===========================================================================
# ROI Generation Script for AFNI Software
#
#---------------------------------------------------------------------------
# This script is designed to generate Region of Interests (ROIs) using the
# AFNI software. It serves as part of the testing suite for assessing the
# functionality and performance of AFNI in creating ROIs. The generated ROIs
# were used to evaluate its accuracy, reproducibility, and robustness, as
# reported in the paper "fMROI: a simple and adaptable toolbox for easy
# region-of-interest creation," currently under review.
#
# To execute the script "genrois_afni.sh," ensure that it and the files
# afni_clustermask.sh, afni_img2mask.sh, afni_cubicmask.sh, and
# afni_spheremask.sh are executable and are in the same directory.
#
# Author: Andre Peres, 2024.
#===========================================================================

./afni_spheremask.sh
./afni_cubicmask.sh
./afni_img2mask.sh
./afni_clustermask.sh
