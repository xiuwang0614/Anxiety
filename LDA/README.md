[Tang, S.Y., et al.,(2020)](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/disorder_subtypes/Tang2020_ASDFactors) provided the specific code, we recommended it for reference.

And this folder recorded the steps taken in our study. 

Using ```addpath(genpath('/path/to/LDA/lib/'))``` to attach relevant files in matlab.

**Before analyses, you need to get the FC of HC and PBD participants.**

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
* The visualziation of Factors was processed by R studio and Python.

#### Step 3: Mediation analyses 
