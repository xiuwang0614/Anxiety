[Tang, S.Y., et al.,(2020)](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/disorder_subtypes/Tang2020_ASDFactors) provided the specific code, we recommended it for reference.

And this folder recorded the steps taken in our study. 

Using ```addpath(genpath('/path/to/LDA/lib/'))``` to attach relevant files in matlab.

**Before analyses, you need to get the FC matrix of HC and PBD participants.**

`FC matrix.mat`: 
The format of file should be N * N * M double in matlab, N is the number of FC nodes, and M is the number of HC gourp or PBD patients. 
The files of HC and PBD group should be separated.

`sub_info_file.csv`: 
The file contains the demograph information of all participants, and the squence should be corresponded to the third dimention of FC matrix.mat

### Usage
#### Step 1: FC2doc
* Folder `step1_FC2doc` contains the functions to convert RSFC data to documents for following steps.
* The z-normalization of PBD patients also completed by this step.

#### Step 2: Estimating latent factors
* the process of estimation should be executed by shell script in the lunix system.
  
#### Step 3: Bootstrapping
* get Bootstrapped Samples.m in matlab
* Bootstrapped Est wrappr.sh in lunix
* computer Bootstrap Zscores.m in matlab
* plot Factors Thresholded_wrapper.m in matlab

**utilities** contained some useful files, which need to be attached.

**polarlda-c-dist** was used to re-compile the source code of polarLDA by the following commands on terminal:
```
# copy to your own code directory to keep the common space clean
cp -aR $/polarlda-c-dist <your_code_dir>
# compile
cd <your_code_dir>/polarlda-c-dist/code; make
```
