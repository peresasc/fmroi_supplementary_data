#bin/bash!

#===========================================================================
# AFNI Spherical ROI Generation
#
# This script generates spherical Region of Interests (ROIs) using the AFNI
# software. It assesses AFNI's functionality and performance in creating
# ROIs, as reported in the paper "fMROI: a simple and adaptable toolbox for
# easy region-of-interest creation," currently under review.
#
# Author: Andre Peres, 2024.
#===========================================================================

srcpath=/home/andre/github/fmroi/etc/test_data/T1.nii.gz
outdir=/media/andre/data8t/fmroi/fmroi_qc/dataset/afni-spheremask

if [ ! -d "$outdir" ]; then
mkdir $outdir
fi

# Number of voxels in each dimension
ni=$(3dinfo -ni $srcpath)
nj=$(3dinfo -nj $srcpath)
nk=$(3dinfo -nk $srcpath)

# Considering that the image is isovoxel, 
# we only take the voxel size od dimension i
vsi=$(3dinfo -adi $srcpath)

sig=("-1" "1")

for rvx in $( seq 100 )
do

r=$(($rvx*${vsi%.*}))

is=${sig[$((RANDOM % 2))]}
it1=$((RANDOM % ($ni/2 - $r - 3)))
it2=$(($ni/2))

i=$(($is*$it1+$it2))

js=${sig[$((RANDOM % 2 + 0))]}
jt1=$((RANDOM % ($nj/2 - $r - 3)))
jt2=$(($nj/2))

j=$(($js*$jt1+$jt2))

ks=${sig[$((RANDOM % 2 + 0))]}
kt1=$((RANDOM % ($nk/2 - $r - 3)))
kt2=$(($nk/2))

k=$(($ks*$kt1+$kt2))


curpos=$(echo $i $j $k)

outpath=$(printf ${outdir}/afni-spheremask_srcimg_t1_radius_%03d_center_x$(($i+1))y$(($j+1))z$(($k+1)).nii $rvx);

# AFNI command for generating a spherical mask
echo $curpos | 3dUndump -srad $r -master $srcpath -prefix $outpath -ijk -

done
