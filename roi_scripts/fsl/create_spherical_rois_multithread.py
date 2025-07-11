#==========================================================================
# FSL Spherical and Cubic ROI Generation
#
#--------------------------------------------------------------------------
# This script generates spherical and Cubic Region of Interests (ROIs) 
# using the FSL package. It assesses FSL's functionality and performance
# in creating ROIs, as reported in the paper "fMROI: a simple and 
# adaptable toolbox for easy region-of-interest creation," currently
# under review.
#
# Author: Igor Vaz, 2024.
#==========================================================================

import os
import random
import argparse
import multiprocessing
from pathlib import Path
import subprocess
from typing import List

import nibabel as nib


def generate_roi_spherical(radius_batch: List[int], t1_image_path: str, output_dir: str) -> None:
    """
    Generate spherical regions of interest (ROIs) of different radii, save them as NIfTI files,
    and store them in the specified output directory.

    Args:
        radius_batch (list): A list of integers specifying the radii of the spherical ROIs to be generated.
        t1_image_path (str): The path to the T1-weighted image.
        output_dir (str): The directory where the ROIs will be saved.

    Returns:
        None
    """
    t1_image = nib.load(t1_image_path)
    t1_image_shape = t1_image.shape
    for radius in radius_batch:
        x, y, z = (
            random.randint(radius, t1_image_shape[0] - radius),
            random.randint(radius, t1_image_shape[1] - radius),
            random.randint(radius, t1_image_shape[2] - radius)
        )
        roi_file = f'roi_radius_{radius}_center_x{x}y{y}z{z}.nii.gz'
        roi_path = Path(output_dir) / roi_file
        command = f'fslmaths {t1_image_path} -mul 0 -add 1 -roi {x} 1 {y} 1 {z} 1 0 1 {roi_path} -odt float && fslmaths {roi_path} -kernel sphere {radius} -fmean -thr 0.0000002 -bin {roi_path} -odt float'
        subprocess.run(command, shell=True, check=True)


def generate_roi_cubic(edge_length_batch: List[int], t1_image_path: str, output_dir: str) -> None:
    """
    Generates cubic regions of interest (ROIs) for a T1-weighted image.

    Args:
        edge_length_batch: A list of integers, each specifying the edge length of a cubic ROI to generate.
        t1_image_path: A string specifying the path to the T1-weighted image to generate ROIs for.
        output_dir: A string specifying the path to the output directory to save the ROIs in.

    Returns:
        None.
    """
    t1_image = nib.load(t1_image_path)
    t1_image_shape = t1_image.shape
    for edge_length in edge_length_batch:
        x, y, z = (
            random.randint(edge_length, t1_image_shape[0] - edge_length),
            random.randint(edge_length, t1_image_shape[1] - edge_length),
            random.randint(edge_length, t1_image_shape[2] - edge_length)
        )
        # https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Fslutils
        roi_file = f'fsl-cubicmask_srcimg_t1_edge_length_{edge_length}_center_x{x}y{y}z{z}.nii.gz'
        roi_path = Path(output_dir) / roi_file
        command = f'fslmaths {t1_image_path} -mul 0 -add 1 -roi {x-edge_length//2} {edge_length} {y-edge_length//2} {edge_length} {z-edge_length//2} {edge_length} 0 1 {roi_path} -odt float'
        subprocess.run(command, shell=True, check=True)


def main():
    parser = argparse.ArgumentParser(description='Generate ROIs')
    parser.add_argument('t1_image_path', help='Path to T1-weighted image')
    parser.add_argument('--output_dir', default='./rois', help='Output directory for ROIs')
    parser.add_argument('--spherical', action='store_true', help='Generate spherical ROIs')
    parser.add_argument('--cubic', action='store_true', help='Generate cubic ROIs')
    parser.add_argument('--batch_size', type=int, default=15, help='Batch size for parallel processing')

    args = parser.parse_args()

    if not args.spherical and not args.cubic:
        parser.error('At least one of --spherical or --cubic is required.')

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    if args.spherical:
        radius = range(1, 101)
        radius_batches = [(radius[i:i+args.batch_size], args.t1_image_path, args.output_dir) for i in range(0, len(radius), args.batch_size)]
        with multiprocessing.Pool(maxtasksperchild=1) as pool:
            pool.starmap(generate_roi_spherical, radius_batches)

    if args.cubic:
        edge_lengths = range(1, 101)
        edge_length_batches = [(edge_lengths[i:i+args.batch_size], args.t1_image_path, args.output_dir) for i in range(0, len(edge_lengths), args.batch_size)]
        with multiprocessing.Pool(maxtasksperchild=1) as pool:
            pool.starmap(generate_roi_cubic, edge_length_batches)


if __name__ == '__main__':
    main()

