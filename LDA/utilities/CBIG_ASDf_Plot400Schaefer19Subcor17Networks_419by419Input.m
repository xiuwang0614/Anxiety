function [corr_mat]=CBIG_ASDf_Plot400Schaefer19Subcor17Networks_419by419Input(corr_mat, scalelim, filename_prefix)
% CBIG_ASDf_Plot400Schaefer19Subcor17Networks_419by419Input(corr_mat, scalelim, filename_prefix)
%
% This function draws the 419x419 correlation matrix (400 cortical ROIs +
% 19 subcortical ROIs). 
% Major networks are separated by thick white grid lines.
% Thin white grid lines separate the breakdowns of major networks and the
% 19 subcortical regions.
% Subcortical structures in the striatum are arranged together in the plot.
% Details on the ordering of the cortical networks and subcortical structures
% are shown below.
%
% This function assumes that the order of the subcortical ROIs in the
% input correlation matrices is in ascending order based on the labels in
% $FREESURFER_HOME/ASegStatsLUT.txt file.
%
% Ordering of major networks and subcortical regions from left to right:
%
%     Major network/Subcortical structures        Sub-network
%     
%     Cortical networks:
%
%     1)  Default:                                TempPar
%                                                 DefaultC 
%                                                 DefaultB
%                                                 DefaultA
%
%     2)  Control:                                ContC
%                                                 ContB 
%                                                 ContA
%
%     3)  Limbic
%
%     4)  SalVentAttn:                            SalVentAttnB
%                                                 SalVentAttnA
%
%     5)  DorsAttn:                               DorsAttnB
%                                                 DorsAttnA
%
%     6)  SomMot:                                 SomMotB
%                                                 SomMotA
%
%     7)  Visual:                                 VisPeri
%                                                 VisCent
%     
%     Subcortical structures:
%     The following 4 striatum structures are arranged together:
%     1)  Accumbens
%     2)  Caudate
%     3)  Pallidum
%     4)  Putamen
%
%     5)  Thalamus
%
%     6)  Amygdala
%     
%     7)  Hippocampus
%
%     8)  Brain Stem
%
%     9)  Diencephalon (Ventral)
%
%    10)  Cerebellum 
%
% Within each sub-network or subcortical region, correlation matrix entries
% start from left hemisphere, then right hemisphere entries (from left to
% right).
%
% Note that the highly "dense" white lines will become more appropriate
% when the figure is saved.
% 
% Input:
%     - corr_mat:
%           419x419 matrix.
%           The order of the subcortical ROIs in this matrix is assumed to 
%           be in ascending order based on the labels in $FREESURFER_HOME/ASegStatsLUT.txt file.
%     - scalelim (optional):
%           Min and max scale limit 
%           If not specified, or specified as [], scale limit is 
%           from -1*max(abs(corr_mat)) to max(abs(corr_mat)).
%     - filename_prefix (optional):
%           Prefix of the saved file name 
%           If not specified, figure will not be saved.
%
% Example:
%       CBIG_ASDf_Plot400Schaefer19Subcor17Networks_419by419Input(corr_mat, [], 'corr_mat')
%       This function will plot corr_mat with max/min scales depending on
%       the maximum absolute value of corr_mat, and save figure as 'corrmat_minsc-1_maxsc1.jpg'.
%
%       CBIG_ASDf_Plot400Schaefer19Subcor17Networks_419by419Input(corr_mat, [-0.5 0.5], [])
%       This function will plot corr_mat with scale limit from -0.5 to 0.5, and will not save figure.
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Check input variables
if exist('scalelim','var')
    if ~isempty(scalelim)
        if size(scalelim,1) > 1
            error('Input argument ''scalelim'' should be a row vector');
        end
    end
end
    
% Get grids
[Index, major_grid, minor_grid, subcor_grid] = LabelsRearrangebyNetwork;
corr_mat = corr_mat(Index,Index);

% Load colormap
load corr_mat_colorscale.mat;


% Plot corr_mat using imagesc
figure; imagesc(corr_mat); set(gcf,'Colormap',rbmap2);
    
% Generate thin grid lines
[xline, yline, ymaj] = generateline(size(corr_mat,1));

xlim(gca,[1 size(corr_mat, 1)]);
ylim(gca,[1 size(corr_mat, 1)]);
postn = get(gca, 'Position');
postn(2) = 0.15;
set(gca, 'Position', postn);

% Set colorbar
hcol=colorbar('peer',gca,'SouthOutside','FontSize',15);
cpos=get(hcol,'Position');
cpos(4)=cpos(4)/4; % Halve the thickness
cpos(3)=cpos(3)*0.75; % Reduce length
cpos(1)=cpos(1) + 0.1; % Move it to the center
cpos(2)=cpos(2) - 0.12; % Move it down outside the plot
set(hcol,'Position',cpos);


% Set color limit 
if ((nargin < 2) || (isempty(scalelim)))
	collim = max(max(abs(corr_mat)));
	scalelim = [-1*collim, 1*collim];
end

set(gca, 'CLim', scalelim);

axis equal;
grid off;
axis([-1 size(corr_mat, 1)+1 -1 size(corr_mat, 1)+1]);
set(gca, 'Visible', 'off')
set(gcf, 'color', 'white');

% Generate major and minor gridlines 
patch(xline(:,subcor_grid), yline(:,subcor_grid),'w', 'edgecolor', 'w', 'Linewidth', 0.005, 'EdgeAlpha', 0.2);
patch(yline(:,subcor_grid), xline(:,subcor_grid),'w', 'edgecolor', 'w', 'Linewidth', 0.005, 'EdgeAlpha', 0.2);
patch(xline(:,minor_grid), yline(:,minor_grid),'w', 'edgecolor', 'w', 'Linewidth', 0.3, 'EdgeAlpha', 0.9);
patch(yline(:,minor_grid), xline(:,minor_grid),'w', 'edgecolor', 'w', 'Linewidth', 0.3, 'EdgeAlpha', 0.9);
patch(xline(:,major_grid), ymaj(:,major_grid),'w', 'edgecolor', 'w', 'Linewidth', 1.1);
patch(ymaj(:,major_grid), xline(:,major_grid),'w', 'edgecolor', 'w', 'Linewidth', 1.1);
    
% save figure
if ((nargin == 3) && ~isempty(filename_prefix))
    filenamefin = [filename_prefix '_minsc' num2str(scalelim(1), '%.1e') '_maxsc' num2str(scalelim(2), '%.1e') '.jpg'];
    set(gcf, 'PaperPositionMode', 'auto','Renderer', 'ZBuffer'); 
    print(gcf, '-djpeg', '-r600', filenamefin);
    close all
end


function [x, y, ymaj] = generateline(n)
x = 1.5:1:n;
x = [ x; x; repmat(nan,1,(n-1)) ];
y = [ 0.5 n+0.5 nan ].';
y = repmat(y,1,(n-1));

ymaj = [ -1 n+1 nan]';
ymaj = repmat(ymaj,1,(n-1));


