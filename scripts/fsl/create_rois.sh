#!/bin/bash

# Path to your python script
SPHERICAL_CUBIC_ROIS_SCRIPT_PATH="/home/igor/proaction/fmroi/create_spherical_rois_multithread.py"

# Call the script with the first set of arguments
python $SPHERICAL_CUBIC_ROIS_SCRIPT_PATH /home/igor/proaction/fmroi/templates/T1.nii.gz --output_dir /home/igor/proaction/fmroi/spherical_rois --spherical

# Call the script with the second set of arguments
python $SPHERICAL_CUBIC_ROIS_SCRIPT_PATH /home/igor/proaction/fmroi/templates/T1.nii.gz --output_dir /home/igor/proaction/fmroi/cubic_rois --cubic
