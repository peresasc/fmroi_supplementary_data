#bin/bash!

#===========================================================================
# AFNI Image Threshold ROI Generation (img2mask)
#
# This script generates Region of Interests (ROIs) based on minimum and
# maximum threshold using the AFNI software. It assesses AFNI's functionality
# and performance in creating ROIs, as reported in the paper "fMROI: a simple
# and adaptable toolbox for easy region-of-interest creation," currently under
# review.
#
# Author: Andre Peres, 2024.
#===========================================================================

srcpath=/home/andre/github/fmroi/etc/test_data/default_mode_association-test_z_FDR_0.01.nii.gz
outdir=/media/andre/data8t/fmroi/fmroi_qc/dataset/afni-img2mask

if [ ! -d "$outdir" ]; then
mkdir $outdir
fi

maxv=$(3dBrickStat -slow $srcpath)

for i in $( seq 100 )
do

at1=$(awk -v m="$maxv" -v n=1 -v seed="$RANDOM" 'BEGIN { srand(seed); for (i=0; i<n; ++i) printf("%.2f\n", rand()*m) }')
t1=$(echo $at1|tr ',' '.')

at2=$(awk -v m="$maxv" -v n=1 -v seed="$RANDOM" 'BEGIN { srand(seed); for (i=0; i<n; ++i) printf("%.2f\n", rand()*m) }')
t2=$(echo $at2|tr ',' '.')

if [ "$(awk -v n1="$t1" -v n2="$t2" 'BEGIN{print (n1>n2)?1:0 }')" -eq 1 ];
then
   aux=$t1
   t1=$t2
   t2=$aux
fi

outpath=$(printf ${outdir}/afni-img2mask_srcimg_dmn_threshold_${t1}_${t2}.nii)

# AFNI command for generating a threshold mask
3dcalc -a $srcpath -expr "within(a,$t1,$t2)" -prefix $outpath

done

# Source image = Complex-shapes
srcpath=/home/andre/github/fmroi/etc/test_data/complex-shapes.nii.gz

tmin=(1.0 3.0 5.0 5.199999 5.399999 5.599999 5.799999 7.0 9.0 11.0)
tmax=(1.0 4.0 5.0 5.20001 5.40001 5.60001 5.80001 8.0 10.0 12.0)

for ii in $( seq ${#tmin[@]} )
do
i=$(($ii-1))
t1=${tmin[$i]}
t2=${tmax[$i]}

outpath=${outdir}/afni-img2mask_srcimg_syndata_threshold_${t1}_${t2}.nii

# AFNI command for generating a threshold mask
3dcalc -a $srcpath -expr "within(a,$t1,$t2)" -prefix $outpath

done
