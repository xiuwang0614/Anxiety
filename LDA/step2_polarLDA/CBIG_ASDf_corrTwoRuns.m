function [corr, avgCorr] = CBIG_ASDf_corrTwoRuns(inputDir_best, inputDir, K, r_best, r)
% [corr, avgCorr] = CBIG_ASDf_corrTwoRuns(inputDir_best, inputDir, K, r_best, r)
%
% This function uses Hungarian matching algorithm to match factors between
% two runs (random initializations), and computes the pairwise correlation between the matched factor
% of the two runs as well as the averaged correlation.
%
% Input:
%     - inputDir_best:
%           Absolute directory to the path where the "reference" solution
%           is saved
%     - inputDir:
%           Absolute directory to the path where the solution to be matched
%           is saved
%     - K:
%           A string indicating the number of factors
%     - r_best:
%           A string. Run number of the "reference" solution
%     - r:
%           A string. Run number of the solution to be matched
%
% Output:
%     - corr:
%           Kx1 vector, where K is the number of factors. Pairwise
%           correlation between the same factor of two runs
%     - avgCorr:
%           Average of corr
%
% Example:
%	[corr, avgCorr] = CBIG_ASDf_corrTwoRuns('~/output','~/output','3','94','10')
%	This function will perform Hungarian matching to match the estimated factors 
%       in run #10 with run #94 for three-factor estimate, and compute the correlation 
%       between the matched factors as well as the averaged correlation.
% 
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Load inputs and compute E(FC patterns | Factor)
% r_best =1;
beta_best = exp(load(fullfile(inputDir_best,sprintf('k%s/r%s', K, num2str(r_best)),'final.beta')));
%beta_best = exp(load(fullfile(inputDir,sprintf('k%d/r%s',3,'1'),'final.beta')));
rho_best = exp(load(fullfile(inputDir_best,sprintf('k%s/r%s', K, num2str(r_best)),'final.rho')));
%rho_best =exp(load(fullfile(inputDir,sprintf('k%d/r%s',3,'1'),'final.rho')));
Mean_best = beta_best.*(2*rho_best-1);
% if find(isnan(Mean_best))
%     %error('Error: Find NaN in corpus.\n')
%     fprintf('Error: Find NaN in corpus.\n');
%     Mean_best(isnan(Mean_best)) = 0;
% end
% if find(isinf(Mean_best))
%     error('Error: Find Inf or -Inf in corpus.\n')
% end
 
beta = exp(load(fullfile(inputDir,sprintf('k%s/r%s', K, num2str(r)),'final.beta')));
% beta = exp(load(fullfile(inputDir,sprintf('k%d/r%s',3, '2'),'final.beta')));
rho = exp(load(fullfile(inputDir,sprintf('k%s/r%s', K, num2str(r)),'final.rho')));
% rho = exp(load(fullfile(inputDir,sprintf('k%d/r%s', 3, '2'),'final.rho')));

Mean = beta.*(2*rho-1);
% if find(isnan(Mean))
%     %error('Error: Find NaN in corpus.\n')
%     fprintf('Error: Find NaN in corpus.\n');
%     Mean(isnan(Mean)) = 0;
% end
% if find(isinf(Mean))
%     error('Error: Find Inf or -Inf in corpus.\n')
% end

%% Reorder factors to obtain the maximal correlation coefficients
K = str2double(K);
% K = str2double(3);
% K=3;
order = CBIG_ASDf_hunMatch(K, Mean_best, Mean);

% Recompute the average correlation with sorted factors
corr = zeros(K, 1);
%corr = zeros(4, 1);
for idx = 1:K  
%for idx = 1:4
	corrMat = corrcoef(Mean_best(idx, :)', Mean(order(idx), :)');
    corr(idx) = corrMat(1, 2); 
end

avgCorr = mean(corr);
