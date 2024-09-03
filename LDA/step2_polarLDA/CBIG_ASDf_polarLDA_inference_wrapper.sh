#!/bin/sh

# Wrapper script to infer factor compositions of new participants with polarLDA model

# Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###################
# Input variables
###################
corpusDir=/media/wang/xiuwang/lda/input/infer/step1_output_infer_dx1.dat # document corpus
outputDir=/media/wang/xiuwang/lda/infer # your output directory
outputName=infFactorComp # your output name
proj_dir=/media/wang/xiuwang/lda
modelDir=${proj_dir}/k4/r13/final # learned model
codeDir=${proj_dir}/lib2/step2_polarLDA# your code directory
infSettings=${code_dir}/CBIG_ASDf_polarLDA_infSettings.txt # inference parameters setting file
#inputDir=$7    # directory of the estimated E(RSFC patterns|Factor) and Pr(Factor|Participant)
factorNum=4 # number of factors

####################################
#                                                                                                                                       
####################################
#ref_run=$(basename $(readlink -f "${inputDir}/k${factorNum}/r*"))
ref_runNum=1
#################
# Run inference
#################
csh ${codeDir}/CBIG_ASDf_polarLDA_inf.csh \
-corpus ${corpusDir} \
-model_dir ${modelDir} \
-factor_num ${factorNum} \
-run_num ${ref_runNum} \
-output_dir ${outputDir} \
-output_name ${outputName} \
-infSettings ${infSettings} \
-code_dir ${codeDir}

########################################
# Write normalized factor compositions
########################################
#gamma_file=${outputDir}/k${factorNum}r${ref_runNum}_${outputName}-gamma.dat
#factorComp_fileName=${outputDir}/${outputName}_k${factorNum}r${ref_runNum}_factorComp.txt
#cd ${codeDir}
#matlab -nodisplay -nosplash -nodesktop -r \
#"clear all;close all;clc; \
#CBIG_ASDf_gamma2table('${gamma_file}','${factorComp_fileName}');exit;"


