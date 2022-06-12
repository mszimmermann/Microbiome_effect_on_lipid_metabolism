function [t,t_amount, t_mean_amount, metNamesMap, gd, useForFitting, curVolumes] = ...
    load_data_from_combined_table(filename, weightfilename)
%% Import data from table file.
dataTable = readtable(filename);
weightTable = readtable(weightfilename);

% combine data into a matrix with time versus species
metNames = dataTable.Properties.VariableNames;
%remove time and group
metNames(cellfun(@(x) isequal(x, 'Time'), metNames))=[];
metNames(cellfun(@(x) isequal(x, 'Group'), metNames))=[];


% create a table for symbiology toolbox
t = sortrows(dataTable, 'Time');

metNamesMap = [metNames' metNames'];
metNamesMap(:,1) = cellfun(@(x) strrep(x, 'FA16_0Mz275', 'D5FA16'),metNamesMap(:,1), 'unif',0);
metNamesMap(:,1) = cellfun(@(x) strrep(x, 'FA16_0Mz301', 'D31FA16'),metNamesMap(:,1), 'unif',0);
% remove colon tissue and iWAT 
%metNamesMap(cellfun(@(x) contains(x, 'iWAT'), metNamesMap(:,1)),:)=[];
metNamesMap(cellfun(@(x) contains(x, '_Colon'), metNamesMap(:,1)),:)=[];

useForFitting = ones(size(metNamesMap,1),1); 
% save curVolumes
curVolumes = ones(size(t,1), length(metNames)); %add extra 0 for time at the beginning

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert concentrations to amounts for fitting
metNames_tissues = metNames;
metNames_tissues = cellfun(@(x) strrep(x, 'FA16_0Mz275_', ''),metNames_tissues, 'unif',0);
metNames_tissues = cellfun(@(x) strrep(x, 'FA16_0Mz301_', ''),metNames_tissues, 'unif',0);
% split group and tissue into two vectors
metNames_group = cellfun(@(x) strsplit(x, '_'), metNames_tissues, 'unif', 0);
metNames_group = cellfun(@(x) x(1), metNames_group, 'unif', 0);
%remove group info
metNames_tissues = cellfun(@(x) replace(x, 'GF_', ''), metNames_tissues, 'unif', 0);
metNames_tissues = cellfun(@(x) replace(x, 'SPF_', ''), metNames_tissues, 'unif', 0);
metNames_tissues = cellfun(@(x) replace(x, 'OMM_', ''), metNames_tissues, 'unif', 0);
% replace OMM with OMM12
metNames_group = cellfun(@(x) strrep(x, 'OMM', 'OMM12'), metNames_group, 'unif', 0);

t_amount = t;
for j=1:length(metNames_tissues)
    weight_i = ismember(weightTable.Tissue, metNames_tissues{j});
    if nnz(weight_i)
        curVol = weightTable{weight_i, metNames_group{j}};
        curVolumes(:,j) = curVol;
        t_amount{:,metNames{j}} = t_amount{:,metNames{j}}*curVol;
    end
   
end
        
t_mean_amount = t_amount;
t_time_unique = unique(t_mean_amount.Time);
for i=1:length(t_time_unique)
    curreplicates = (t_mean_amount.Time == t_time_unique(i));
    for j=1:length(metNames)
        t_mean_amount{curreplicates,metNames{j}} = nanmean(t_amount{curreplicates,metNames{j}});
    end
end
t_mean_amount = unique(t_mean_amount);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert data to SimBiology format
gd = groupedData(t_amount);
useForFitting = find(useForFitting);
