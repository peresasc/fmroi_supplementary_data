#bin/bash!

#===========================================================================
# AFNI Clustering ROI Generation
#
# This script generates Region of Interests (ROIs) using the 3dClusterize
# algorithm in AFNI. It assesses AFNI's functionality and performance in
# creating ROIs, as reported in the paper "fMROI: a simple and adaptable
# toolbox for easy region-of-interest creation," currently under review.
#
# Author: Andre Peres, 2024.
#===========================================================================

srcpath=/home/andre/github/fmroi/etc/test_data/64-spheres.nii.gz
outdir=/media/andre/data8t/fmroi/fmroi_qc/dataset/afni-clustermask_map

if [ ! -d "$outdir" ]; then
mkdir $outdir
fi

maxv=$(3dBrickStat -slow $srcpath)

tmin=(0.1 17.0 33.0)
tmax=($maxv 32 $maxv)

mincsz=(1 33 123)

for ii in $( seq ${#tmin[@]} )
do
i=$(($ii-1))
t1=${tmin[$i]}
t2=${tmax[$i]}
cs=${mincsz[$i]}

outpath=${outdir}/afni-clustermask_srcimg_spheres_threshold_${t1}_${t2}_mincsz_${cs}.nii

# AFNI command for generating the cluster masks
3dClusterize                   \
    -inset $srcpath           \
    -ithr 0                    \
    -idat 0                    \
    -NN 1                      \
    -within_range $t1 $t2     \
    -clust_nvox $cs            \
    -no_1Dformat           \
    -pref_map $outpath      \
    
done
