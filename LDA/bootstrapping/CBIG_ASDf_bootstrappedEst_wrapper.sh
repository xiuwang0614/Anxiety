#!/bin/bash
# 
# Wrapper script to run polarLDA estimate on bootstrapped samples
# 
# Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Take input variables
input_dir=/media/wang/OpenData/4LDA/output/bootstrap/data # absolute directory where the docs are saved 
out_dir=/media/wang/OpenData/4LDA/output/bootstrap/output
K=4  # number of factors
N=100  # number of resamples
#cluster=$5 # cluster name

proj_dir=/media/wang/OpenData/4LDA
code_dir=${proj_dir}/lib2/step2_polarLDA
model_name=${proj_dir}/output/visualization/k4/r13/final/final
#run_files = ${code_dir}/run_files_1.txt
inf_settings=${code_dir}/CBIG_ASDf_polarLDA_infSettings.txt

for (( i=1; i<=${N}; i++ ))
do
    corpusDir_est=${input_dir}/resampled_${i}/dx1.dat
    output_dir=${out_dir}/resampled_${i}
    mkdir -p ${output_dir}
    output_dir_step2a=${output_dir}/estimate
    #mkdir -p ${output_dir_step2a}

    sh ${code_dir}/CBIG_ASDf_polarLDA_est_initFromModel.sh \
    -d ${corpusDir_est} \
    -t ${inf_settings} \
    -k ${K} \
    -i model \
    -m ${model_name} \
    -p ${code_dir} \
    -o ${output_dir_step2a} 
#    -q ${cluster}

done
