% Parse gene expression data and merge in one table
dataFolder = '.\Data\expression_atlas\Affymetrix\';
%dataFolder = '.\Data\expression_atlas\RNAseq\';
dataSubfolders = dir(dataFolder);
dataSubfolder_names = cell(length(dataSubfolders),1);
for i=1:length(dataSubfolders)
    if length(dataSubfolders(i).name)>2
        dataSubfolder_names{i} = dataSubfolders(i).name;
    end
end
dataSubfolder_names(cellfun(@(x) isempty(x), dataSubfolder_names))=[];


totalAssays = cell(10000,2);
totalContrasts = cell(10000,2);
idx = 1;
cidx = 1;
totalFCtable = [];
totalPtable = [];
% read files from each folder and get the fold changes and the experiment
% info
for i=1:length(dataSubfolder_names)
    try

        curFolder = dir([dataFolder, dataSubfolder_names{i}]);
        file_analytics = cellfun(@(x) contains(x, 'analytics'), {curFolder(:).name});
        file_analytics = curFolder(file_analytics).name;
        % get description
        file_desc = cellfun(@(x) contains(x, 'configuration'), {curFolder(:).name});
        file_desc = curFolder(file_desc).name;
        
        % read fold changes
        curDataTable = readtable([dataFolder, dataSubfolder_names{i}, filesep, ...
            file_analytics], 'FileType', 'text');
        col_gene = cellfun(@(x) contains(x, 'GeneID'), curDataTable.Properties.VariableNames);
        col_pval = cellfun(@(x) contains(x, 'p_value'), curDataTable.Properties.VariableNames);
        col_log2fc = cellfun(@(x) contains(x, 'log2foldchange'), curDataTable.Properties.VariableNames);
        col_names = curDataTable.Properties.VariableNames;
        curFCtable = curDataTable(:, col_names(col_gene | col_log2fc));
        % add experiment ID to columns
        edit_cols =  cellfun(@(x) contains(x, 'log2foldchange'), curFCtable.Properties.VariableNames);      
        curFCtable.Properties.VariableNames(edit_cols) = cellfun(@(x) strcat(dataSubfolder_names{i}, '_', x),...
            curFCtable.Properties.VariableNames(edit_cols), 'unif', 0);
        % keep only unique GeneIDs
        [~, uniqueidx] = unique(curFCtable.GeneID);
        curFCtable = curFCtable(uniqueidx,:);
        curFCtable.Properties.RowNames = curFCtable.GeneID;
        curFCtable.GeneID = [];
        % get pvalue table
        curPtable = curDataTable(:, col_names(col_gene | col_pval));
        % add experiment ID to columns
        edit_cols =  cellfun(@(x) contains(x, 'p_value'), curPtable.Properties.VariableNames);      
        curPtable.Properties.VariableNames(edit_cols) = cellfun(@(x) strcat(dataSubfolder_names{i}, '_', x),...
            curPtable.Properties.VariableNames(edit_cols), 'unif', 0);
        % keep only unique GeneIDs
        [~, uniqueidx] = unique(curPtable.GeneID);
        curPtable = curPtable(uniqueidx,:);
        curPtable.Properties.RowNames = curPtable.GeneID;
        curPtable.GeneID = [];

        if isempty(totalFCtable)
            totalFCtable = curFCtable;
            totalPtable = curPtable;
        else
            testsize = height(totalFCtable);
           
            totalFCtable = outerjoin(totalFCtable,curFCtable,'Keys','Row');
            totalPtable = outerjoin(totalPtable,curPtable,'Keys','Row');
            testsizemerged = height(totalFCtable);
            if (testsizemerged - testsize) > height(curFCtable)
                break
            end
            
        end
        
     % get configuration information
        curConfig = readstruct([dataFolder, dataSubfolder_names{i}, filesep, ...
            file_desc]);

        try
            % get list of assays
            curAssays = curConfig.analytics.assay_groups.assay_group;
            curAssays_id = {curAssays(:).idAttribute};
            curAssays_label = {curAssays(:).labelAttribute};
            curAssays_id = cellfun(@(x) strcat(dataSubfolder_names{i}, '_', x), curAssays_id, 'unif', 0);
            %add to total
            totalAssays(idx:idx+length(curAssays_id)-1,:) = [curAssays_id' curAssays_label'];
            idx = idx+length(curAssays_id);
        catch
            disp(['Error reading ', dataSubfolder_names{i}, ' array info\n'])
        end
        %%%%%%%%%%%%%%%%%%%%%%
        % get list of contrasts
        try
            curContrasts = curConfig.analytics.contrasts.contrast;
            curContrast_id = {curContrasts(:).idAttribute};
            curContrast_label = {curContrasts(:).name};
            curContrast_id = cellfun(@(x) strcat(dataSubfolder_names{i}, '_', x), curContrast_id, 'unif', 0);
            %add to total
            totalContrasts(cidx:cidx+length(curContrast_id)-1,:) = [curContrast_id' curContrast_label'];
            cidx = cidx+length(curContrast_id);
        catch
            disp(['Error reading ', dataSubfolder_names{i}, ' contrasts info\n'])
        end

    catch
        disp(['Error reading ', dataSubfolder_names{i}, ' \n'])
    end
end
totalContrasts(cidx:end,:) = [];
totalAssays(idx:end,:) = [];


totalContrasts_only = cell(10000,2);
cidx = 1;
% read only contrasts as some were not read due to errors in assay reading
for i=1:length(dataSubfolder_names)
    try

        curFolder = dir([dataFolder, dataSubfolder_names{i}]);
        % get description
        file_desc = cellfun(@(x) contains(x, 'configuration'), {curFolder(:).name});
        file_desc = curFolder(file_desc).name;
        
        %%%%%%%%%%%%%%%%%%%%%%
        % get list of contrasts
        curConfig = readstruct([dataFolder, dataSubfolder_names{i}, filesep, ...
            file_desc]);

            curContrasts = curConfig.analytics.contrasts.contrast;
            curContrast_id = {curContrasts(:).idAttribute};
            curContrast_label = {curContrasts(:).name};
            curContrast_id = cellfun(@(x) strcat(dataSubfolder_names{i}, '_', x), curContrast_id, 'unif', 0);
            %add to total
            totalContrasts_only(cidx:cidx+length(curContrast_id)-1,:) = [curContrast_id' curContrast_label'];
            cidx = cidx+length(curContrast_id);
    catch
        disp(['Error reading ', dataSubfolder_names{i}, ' \n'])
    end
end
totalContrasts_only(cidx:end,:) = [];
% clean up the table as many genes seem to be appearing multiple times
% totalFCtable_cleaned = array2table(nan(length(unique(totalFCtable.GeneID)), ...
%                                     length(totalFCtable.Properties.VariableNames)-2),...
%                             'VariableNames', totalFCtable.Properties.VariableNames(3:end),...
%                             'RowNames', unique(totalFCtable.GeneID));
% 
% for i=3:length(totalFCtable.Properties.VariableNames)
%     curData = totalFCtable{:, totalFCtable.Properties.VariableNames(i)};
%     curGenes = totalFCtable{:,1};
%     curnan = isnan(curData);
%     % remove nans
%     curData(curnan) = [];
%     curGenes(curnan) = [];
%     % get the correct indices in the cleaned table
%     [~, ~, curidx] = intersect(curGenes, totalFCtable_cleaned.Properties.RowNames, 'stable');
%     totalFCtable_cleaned{curidx, i-2} = curData;
% end
% 
% % save table to file
% writetable(totalFCtable_cleaned, '.\Output\totalFCtable_expression_atlas.csv',...
%     'WriteRowNames', 1);
% 
% % clean up the pval table as many genes seem to be appearing multiple times
% totalPtable_cleaned = array2table(nan(length(unique(totalPtable.GeneID)), ...
%                                     length(totalPtable.Properties.VariableNames)-2),...
%                             'VariableNames', totalPtable.Properties.VariableNames(3:end),...
%                             'RowNames', unique(totalPtable.GeneID));
% 
% for i=3:length(totalPtable.Properties.VariableNames)
%     curData = totalPtable{:, totalPtable.Properties.VariableNames(i)};
%     curGenes = totalPtable{:,1};
%     curnan = isnan(curData);
%     % remove nans
%     curData(curnan) = [];
%     curGenes(curnan) = [];
%     % get the correct indices in the cleaned table
%     [~, ~, curidx] = intersect(curGenes, totalPtable_cleaned.Properties.RowNames, 'stable');
%     totalPtable_cleaned{curidx, i-2} = curData;
% end
% 
% % save table to file
% writetable(totalPtable_cleaned, '.\Output\expression_atlas\totalPtable_expression_atlas.csv',...
%     'WriteRowNames', 1);

totalFCtable.Properties.RowNames = totalFCtable.GeneID;
totalFCtable.GeneID = [];
writetable(totalFCtable, '.\Output\expression_atlas\totalFCtable_expression_atlas_affymetrix.csv',...
    'WriteRowNames', 1);

% writetable(totalFCtable, '.\Output\expression_atlas\totalFCtable_expression_atlas.csv',...
%     'WriteRowNames', 1);
totalPtable.Properties.RowNames = totalPtable.GeneID;
totalPtable.GeneID = [];
writetable(totalPtable, '.\Output\expression_atlas\totalPtable_expression_atlas_affymetrix.csv',...
    'WriteRowNames', 1);

% writetable(totalPtable, '.\Output\expression_atlas\totalPtable_expression_atlas.csv',...
%     'WriteRowNames', 1);

% write contrasts and assay tables
totalContrasts_table = cell2table(totalContrasts_only, 'VariableNames', {'ContrastID', 'ContrastName'});
totalAssays_table = cell2table(totalAssays, 'VariableNames', {'AssayID', 'AssayName'});
writetable(totalContrasts_table, '.\Output\expression_atlas\totalContrasts_expression_atlas_affymetrix.csv');
writetable(totalAssays_table, '.\Output\expression_atlas\totalAssays_expression_atlas_affymetrix.csv');
% writetable(totalContrasts_table, '.\Output\expression_atlas\totalContrasts_expression_atlas.csv');
% writetable(totalAssays_table, '.\Output\expression_atlas\totalAssays_expression_atlas.csv');

% save table with ensemble gene id to gene name mapping
[~, ~, curidx] = intersect(totalPtable_cleaned.Properties.Rownames, totalPtable.GeneID, 'stable');

totalGeneID_Names = totalPtable{curidx, 1:2};
totalGeneID_Names = cell2table(totalGeneID_Names, 'VariableNames', {'GeneID', 'GeneName'});
writetable(totalGeneID_Names, '.\Output\expression_atlas\totalGeneID_Names.csv');



