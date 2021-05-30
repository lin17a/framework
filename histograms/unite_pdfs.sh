#!/bin/bash

custom_string='original'

histograms="$(ls '/nfs/dust/cms/user/meyerlin/ba/framework/histograms/m400_w20_pseudo_scalar_res_'"$custom_string"'_generated_vs_reweighted')"

echo $histograms

#for hist in $() 

nhists=`echo "scale=0; ${#MASSES[@]} - 1" | bc`

echo $nhists
#for CONFIG in `seq 0 ${NCONFIGS}`; do
