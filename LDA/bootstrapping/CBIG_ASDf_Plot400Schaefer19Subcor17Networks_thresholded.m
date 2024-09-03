function [corr_mat_parcelOrder,corr_mat] = CBIG_ASDf_Plot400Schaefer19Subcor17Networks_thresholded(corr_mat_row, ...
 scalelim, filename_prefix, blk_indices)
% corr_mat_parcelOrder = CBIG_ASDf_Plot400Schaefer19Subcor17Networks_thresholded(corr_mat_row, 
% scalelim, filename_prefix, blk_indices)
%
% Plots thresholded 419x419 matrix masked by input blk_indices.
% For details about the 419x419 matrix, refer to
% CBIG_ASDf_Plot400Schaefer19Subcor17Networks.
% 
% Input:
%     - corr_mat_row:
%           1 x 87571 row vector, lower triangular entries (exc. diagonal) in 419x419
%           correlation matrix.           
%     - scalelim (optional):
%           Min and max scale limit 
%           If not specified, or specified as [], scale limit is 
%           from -1*max(abs(corr_mat)) to max(abs(corr_mat)).
%     - filename_prefix (optional):
%           Prefix of the saved file name 
%           If not specified, figure will not be saved.
%     - blk_indices (optional):
%           1 x 171 (18x19/2) row vector of boolean values, indicating
%           whether the values in lower triangular entries (inc. diagonal)
%           in 18x18 within/between network matrix are significant or not.
%
% Output:
%     - corr_mat_parcelOrder:
%           1 x 87571 row vector, masked input corr_mat_row based on
%           blk_indices. i.e., if blk_indices(i,j) = true,
%           corr_mat_parcelOrder(i,j) = 0, else corr_mat_parcelOrder(i,j) =
%           corr_mat_row(i,j).
%
% Example:
%       CBIG_ASDf_Plot400Schaefer19Subcor17Networks_thresholded(corr_mat_row, 
%           [-1.6e-5 1.6e-5], 'factor1_thresholded', blk_indices)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% From row vector to 419x419 matrix
lowerTri = tril(ones(416),-1);
lowerTri(lowerTri ~= 0) = corr_mat_row;
corr_mat = lowerTri + lowerTri';

% Convert blk_indices to 18x18 matrix
lowerTr = tril(ones(29));
lowerTr(lowerTr ~= 0) = blk_indices;
blk_indices_mat = lowerTr + lowerTr';
blk_indices_mat = logical(blk_indices_mat);

% Get grids
[Index, major_grid, minor_grid, subcor_grid, Index_back] = LabelsRearrangebyNetwork;
corr_mat = corr_mat(Index,Index);

blk_grid = [unique(minor_grid) 416];
num_blks = length(blk_grid);

% Mask values within/between blocks
for j = 1:num_blks
    if j == 1
        y_start = 1;
    else
        y_start = blk_grid(j-1) + 1;
    end
    y_end = blk_grid(j);
    
    for i = 1:num_blks
        if i == 1
            x_start = 1;
        else
            x_start = blk_grid(i-1) + 1;
        end
        x_end = blk_grid(i);
        
        corr_mat(x_start:x_end,y_start:y_end) = corr_mat(x_start:x_end,y_start:y_end) * blk_indices_mat(i,j);
    end
end

% Back to original ordering
corr_mat_parcelOrder = corr_mat(Index_back, Index_back);

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
axis([-5 size(corr_mat, 1)+5.5 -5 size(corr_mat, 1)+5.5]);
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
if ((nargin > 2) && ~isempty(filename_prefix))
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

ymaj = [ -5 n+5.5 nan]';
ymaj = repmat(ymaj,1,(n-1));



function [Index, major_grid, minor_grid, subcor_grid, Index_back] = LabelsRearrangebyNetwork

%load original cortical networks labels
% networks_path = fullfile('I:\lda\mask\');
% lh_annot_file = [networks_path 'lh.Schaefer2018_400Parcels_7Networks_order.annot'];
% rh_annot_file = [networks_path 'rh.Schaefer2018_400Parcels_7Networks_order.annot'];
% [lh_vertex_labels, lh_colortable] = CBIG_read_annotation(lh_annot_file);
% [rh_vertex_labels, rh_colortable] = CBIG_read_annotation(rh_annot_file);
% 
% lh_label = lh_colortable.struct_names(2:end);
% rh_label = rh_colortable.struct_names(2:end);
% for i = 1:numel(rh_label)
%     if contains(rh_label{i}, '_LH_')
%         rh_label{i} = strrep(rh_label{i}, '_LH_', '_RH_');
%     end
% end
% % hard-coded, assuming the original input subcortical labels are in ascending order
% subcor_labelname =  {'Hippocampus_Left';'Hippocampus_Right';'HippocampusPara_Left';'HippocampusPara_Right';...
%     'Amygdala_Left';'Amygdala_Right';'Caudate_Left';'Caudate_Right';...
%     'Putamen_Left';'Putamen_Right';'Pallidum_Left';'Pallidum_Right';...
%     'Accumbens_Left';'Accumbens_Right';'Thalamus_Left';'Thalamus_Right'};
% 
% major_grid = [];
% minor_grid = [];
% subcor_grid = [];
% %major_acc_index = [1, 4, 7, 8, 10, 12, 14, 16];
% major_acc_index = [3,5,7,8,11,13,15,17];
% subcor_acc_index = 18:25;
% 
% lhrh_label = {lh_label{:} rh_label{:}}';
% %lhrh_label = load('E:\LDA\mask\lhrh_lable.mat');
% %lhrh_label = lhrh_label.lhrh_label;
% all_label = [lhrh_label; subcor_labelname];
% all_netw = cell(numel(all_label),1);
% all_subnet = cell(numel(all_label),1);
T = readtable('I:\lda\mask\Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.Centroid_RAS.csv');
T= T.ROIName;
% hard-coded, assuming the original input subcortical labels are in ascending order
subcor_labelname =  {'Hippocampus_Left';'Hippocampus_Right';'HippocampusPara_Left';'HippocampusPara_Right';...
    'Amygdala_Left';'Amygdala_Right';'Caudate_Left';'Caudate_Right';...
    'Putamen_Left';'Putamen_Right';'Pallidum_Left';'Pallidum_Right';...
    'Accumbens_Left';'Accumbens_Right';'Thalamus_Left';'Thalamus_Right'};

major_grid = [];
minor_grid = [];
subcor_grid = [];
%major_acc_index = [1, 4, 7, 8, 10, 12, 14, 16];
%f2-4:
% major_acc_index = [8,15,17,23,26,27,28];
%f1
major_acc_index = [1,9,16,19,20,26,28];
subcor_acc_index = 29:36;
lhrh_label = T;
% lhrh_label = {lh_label{:} rh_label{:}}';
%lhrh_label = load('E:\LDA\mask\lhrh_lable.mat');
%lhrh_label = lhrh_label.lhrh_label;
all_label = [lhrh_label; subcor_labelname];
all_netw = cell(numel(all_label),1);
all_subnet = cell(numel(all_label),1);

for i = 1:numel(all_label)
    % cortical
    if i <= numel(lhrh_label)
        netw = textscan(char(all_label(i,1)), '%s %s %s %s', 'delimiter', '_');
        tmp = netw{1,3};
        all_netw(i,1) = tmp;
        subnet = netw{1,4};
        if ~isempty(subnet)
            all_subnet(i,1) = subnet;
        end
        
    % subcortical    
    else
        netw = textscan(char(all_label(i,1)),'%s %s', 'delimiter', '_');
        all_netw(i,1) = netw{1,1};       
    end
end

% arrange new labels based on template order
% tmplate = {'TempPar'; 'DefaultC'; 'DefaultB';'DefaultA'; 'ContC'; 'ContB'; 'ContA'; ...
%     'SalVentAttnB'; 'SalVentAttnA'; 'DorsAttnB'; 'DorsAttnA'; 'SomMotB'; 'SomMotA'; 'VisPeri'; 'VisCent'; ...
%     'LimbicA'; 'LimbicB';...
%     'Hippocampus';'HippocampusPara';...
%     'Amygdala';'Caudate';'Putamen';'Pallidum';...
%     'Accumbens';'Thalamus'};
% tmplate = {'ContA';  'ContC'; 'ContB';...
%     'SalVentAttnB'; 'SalVentAttnA';'LimbicB';'LimbicA'; ...
%     'TempPar';'DefaultC'; 'DefaultA';'DefaultB';
%     'VisPeri'; 'VisCent';'SomMotB'; 'SomMotA'; ...
%     'DorsAttnA';'DorsAttnB'; ...
%     'Hippocampus';'HippocampusPara';...
%     'Amygdala';'Caudate';'Putamen';'Pallidum';'Accumbens';'Thalamus'};
% tmplate = {'ContA';'ContC'; 'ContB'; 'SalVentAttnB';'SalVentAttnA';...
%     'LimbicB';'LimbicA';...    
%     'TempPar';'DefaultA';'DefaultC';'DefaultB'; ...
%     'DorsAttnB';'DorsAttnA'; ...
%     'VisPeri'; 'VisCent'; 'SomMotB';'SomMotA'; ...
%     'Hippocampus';'HippocampusPara';...
%     'Amygdala';'Caudate';'Putamen';'Pallidum';'Accumbens';'Thalamus'};
%tmplate2 = {'TempPole'; 'OFC'};
% tmplate = {'ContPar';'ContTemp';'ContOFC';'ContPFCl';'ContPFCv';'ContpCun';'ContCing';'ContPFCmp';...
% 'SalVentAttnParOper';'SalVentAttnTempOcc';'SalVentAttnTempOccPar';'SalVentAttnFrOperIns';'SalVentAttnPFCl';'SalVentAttnMed';'SalVentAttnPrC';...
% 'LimbicOFC';'LimbicTempPole';...
% 'DefaultTemp';'DefaultPar';'DefaultPFC';'DefaultPFCv';'DefaultPFCdPFCm';'DefaultpCunPCC';...
% 'DorsAttnPost';'DorsAttnFEF';'DorsAttnPrCv';...
% 'Vis';'SomMot';...
% 'Hippocampus';'HippocampusPara';...
% 'Amygdala';'Caudate';'Putamen';'Pallidum';'Accumbens';'Thalamus'};
tmplate = {'SomMot';'ContPar';'ContTemp';'ContOFC';'ContPFCl';'ContPFCv';'ContpCun';'ContCing';'ContPFCmp';...
'SalVentAttnParOper';'SalVentAttnTempOcc';'SalVentAttnTempOccPar';'SalVentAttnFrOperIns';'SalVentAttnPFCl';'SalVentAttnMed';'SalVentAttnPrC';...
'DorsAttnPost';'DorsAttnFEF';'DorsAttnPrCv';...
'Vis';...
'DefaultTemp';'DefaultPar';'DefaultPFC';'DefaultPFCv';'DefaultPFCdPFCm';'DefaultpCunPCC';...
'LimbicOFC';'LimbicTempPole';...
'Hippocampus';'HippocampusPara';...
'Amygdala';'Caudate';'Putamen';'Pallidum';'Accumbens';'Thalamus'};
newlabel = [];

curInd = 0;
for j = 1:numel(tmplate)
    ind = find(strcmp(all_netw,tmplate(j)));
    
%    % for Limbic networks, further separate the networks into TempPole and OFC
%    if sum(strcmp('Limbic', tmplate(j))) ~= 0
%        ind2 = [];
%        for s = 1:numel(tmplate2)
%            if sum(strcmp(tmplate2(s), all_subnet(ind))) ~= 0
%                ind2 = [ind2; ind(find(strcmp(tmplate2(s), all_subnet(ind))))];
%                minor_grid = [minor_grid curInd+size(ind2,1)];
%            end
%            
%        end
%        if numel(ind) ~= numel(ind2)
%            disp('Wrong Index')
%        end
%        ind = ind2;
%    end

    curInd = curInd+size(ind,1);
    if (j~=numel(tmplate))
        if (any(j==subcor_acc_index))
            subcor_grid = [subcor_grid curInd];
        else
            minor_grid = [minor_grid curInd];
        end
    end
    if (any(j==major_acc_index))
        major_grid = [major_grid curInd];
    end
    
    newlabel = cat(1, newlabel,all_label(ind));
end

% create indexing from old labeling to new labeling
Index = zeros(numel(all_label),1);
Index_back = zeros(numel(all_label),1);
for k = 1:numel(all_label)
    
    Index(k) = find(strcmp(newlabel(k),all_label));
    Index_back(Index(k)) = k;
end

fin_label = all_label(Index);
