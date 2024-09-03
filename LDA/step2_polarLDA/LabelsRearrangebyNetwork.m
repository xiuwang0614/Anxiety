function [Index, major_grid, minor_grid, subcor_grid] = LabelsRearrangebyNetwork

% %load original cortical networks labels
% networks_path = fullfile('I:\lda\mask\0518');
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
major_acc_index = [8,15,17,23,26,27,28];
%f1
%major_acc_index = [1,9,16,19,20,26,28];
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
% tmplate = { 'VisPeri'; 'VisCent';'DefaultC'; 'DefaultB';'DefaultA'; ...
%     'TempPar';'DorsAttnB'; 'DorsAttnA';'SomMotB'; 'SomMotA';...
%      'LimbicA'; 'LimbicB';'Hippocampus';'HippocampusPara';...
%     'Amygdala';'Caudate';'Putamen';'Pallidum';...
%     'Accumbens';'Thalamus';...
%     'ContC'; 'ContB'; 'ContA';'SalVentAttnB'; 'SalVentAttnA'};
% tmplate = {'ContA';'ContC'; 'ContB'; 'SalVentAttnB';'SalVentAttnA';...
%     'LimbicB';'LimbicA';...    
%     'TempPar';'DefaultA';'DefaultC';'DefaultB'; ...
%     'DorsAttnB';'DorsAttnA'; ...
%     'VisPeri'; 'VisCent'; 'SomMotB';'SomMotA'; ...
%     'Hippocampus';'HippocampusPara';...
%     'Amygdala';'Caudate';'Putamen';'Pallidum';'Accumbens';'Thalamus'};
%tmplate2 = {'TempPole'; 'OFC'};
% F2-4
tmplate = {'ContPar';'ContTemp';'ContOFC';'ContPFCl';'ContPFCv';'ContpCun';'ContCing';'ContPFCmp';...
'SalVentAttnParOper';'SalVentAttnTempOcc';'SalVentAttnTempOccPar';'SalVentAttnFrOperIns';'SalVentAttnPFCl';'SalVentAttnMed';'SalVentAttnPrC';...
'LimbicOFC';'LimbicTempPole';...
'DefaultTemp';'DefaultPar';'DefaultPFC';'DefaultPFCv';'DefaultPFCdPFCm';'DefaultpCunPCC';...
'DorsAttnPost';'DorsAttnFEF';'DorsAttnPrCv';...
'Vis';'SomMot';...
'Hippocampus';'HippocampusPara';...
'Amygdala';'Caudate';'Putamen';'Pallidum';'Accumbens';'Thalamus'};
%F1
% tmplate = {'SomMot';'ContPar';'ContTemp';'ContOFC';'ContPFCl';'ContPFCv';'ContpCun';'ContCing';'ContPFCmp';...
% 'SalVentAttnParOper';'SalVentAttnTempOcc';'SalVentAttnTempOccPar';'SalVentAttnFrOperIns';'SalVentAttnPFCl';'SalVentAttnMed';'SalVentAttnPrC';...
% 'DorsAttnPost';'DorsAttnFEF';'DorsAttnPrCv';...
% 'Vis';...
% 'DefaultTemp';'DefaultPar';'DefaultPFC';'DefaultPFCv';'DefaultPFCdPFCm';'DefaultpCunPCC';...
% 'LimbicOFC';'LimbicTempPole';...
% 'Hippocampus';'HippocampusPara';...
% 'Amygdala';'Caudate';'Putamen';'Pallidum';'Accumbens';'Thalamus'};


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

for k = 1:numel(all_label)
    
    Index(k) = find(strcmp(newlabel(k),all_label));
end

fin_label = all_label(Index);
