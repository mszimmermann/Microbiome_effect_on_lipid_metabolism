% analyze proteomics data from the liver and try to find conditions from
% the public databases (Expression Atlas) to identify under which
% conditions these genes are also differentially expressed

% load proteomics data
fileName = '.\Data\proteome_data\Proteomics_GF_OMM12-SPF_Liver_Data only_221024.xlsx';
figFolder = '.\Figures\';
outputFolder = '.\Output\';

protTable = readtable(fileName);
        
% get conditions from table columns
tableColumns = protTable.Properties.VariableNames;
tableColumns_groups = cellfun(@(x) (contains(x, '_') &...
                                                 ~contains(lower(x), 'mean') &...
                                                 ~contains(lower(x), 'vs')), tableColumns);
tableColumns_groups_names = cellfun(@(x) x(1:strfind(x,'_')-1), tableColumns, 'unif', 0);
tableColumns_groups_unique = unique(tableColumns_groups_names);

% perform differential analysis of proteins between colonized and GF groups
fcMatrix_OMM_SPF_vsGF = zeros(height(protTable), 2);
pMatrix_OMM_SPF_vsGF = zeros(height(protTable), 2);
fdrMatrix_OMM_SPF_vsGF = zeros(height(protTable), 2);

contrast_groups = {'OMM12', 'SPF'};
ref_group = 'GF';
for i=1:length(contrast_groups)
    curFC = nanmean(protTable{:, ismember(tableColumns_groups_names,...
                                  contrast_groups{i}) &...
                                  tableColumns_groups},2)-...
            nanmean(protTable{:, ismember(tableColumns_groups_names,...
                                  ref_group) &...
                                  tableColumns_groups},2);

    [~, pMatrix] = ttest2(protTable{:, ismember(tableColumns_groups_names,...
                                  contrast_groups{i}) &...
                                  tableColumns_groups},...
                          protTable{:, ismember(tableColumns_groups_names,...
                                  ref_group) &...
                                  tableColumns_groups},...
                          'Dim', 2, 'VarType', 'equal');

     pFDR = mafdr(pMatrix, 'bhfdr', 1);

     fcMatrix_OMM_SPF_vsGF(:,i) = curFC;
     pMatrix_OMM_SPF_vsGF(:,i) = pMatrix;
     fdrMatrix_OMM_SPF_vsGF(:,i) = pFDR;
end

% select proteins based on thresholds in both groups
pFDRthreshold = 0.05;
fcThreshold = log2(1.5);
select_proteins = ((abs(fcMatrix_OMM_SPF_vsGF(:,1))>=fcThreshold) &...
                   (fdrMatrix_OMM_SPF_vsGF(:,1)<=fcThreshold) &...
                   (abs(fcMatrix_OMM_SPF_vsGF(:,2))>=fcThreshold) &...
                   (fdrMatrix_OMM_SPF_vsGF(:,2)<=fcThreshold));
                   
% get IDs of the selected proteins
select_proteins_geneids = protTable{select_proteins, 12};
select_proteins_uniprotids = protTable{select_proteins, 'Accession'};

proteins_uniprotids = protTable{:, 'Accession'};

% convert protein IDs from uniprot to ensembl
%load uniprot conversion table
uniprot_file = '.\Data\UNIPROT\MOUSE_10090_idmapping.dat';
uniprot_file = 'Z:\mazimmer\Projects\MZK006_Collaborations\MZK006D_Microbiome_effect_on_lipid_metabolism\mzk006d_microbiome_effect_on_lipid_metabolism\Data\UNIPROT\MOUSE_10090_idmapping.dat';
uniprot_table = readtable(uniprot_file, 'delim', '\t');

%convert one by one all IDs to ENSEMBLEGenes

conversion_table = (cellfun(@(x) ismember(x, select_proteins_uniprotids), uniprot_table{:,1}) &...
        cellfun(@(x) ismember(x, {'Ensembl'}), uniprot_table{:,2}));

select_proteins_ensembleids = uniprot_table(conversion_table,:);


conversion_table = (cellfun(@(x) ismember(x, proteins_uniprotids), uniprot_table{:,1}) &...
        cellfun(@(x) ismember(x, {'Ensembl'}), uniprot_table{:,2}));

proteins_ensembleg = uniprot_table(conversion_table,:);
% leave only unique ensemble genes
[~, uniqueidx] = unique(proteins_ensembleg{:,3});
proteins_ensembleg = proteins_ensembleg{uniqueidx,:};

% write fold changes and p-values for ensemble converted proteins
[~, idxFC, idxG] = intersect(proteins_uniprotids, proteins_ensembleg(:,1));
fcTable_ensembleg = fcMatrix_OMM_SPF_vsGF(idxFC,:);
fcTable_ensembleg = array2table(fcTable_ensembleg, 'RowNames', proteins_ensembleg(idxG,3),...
    'VariableNames', strcat('log2FC_',contrast_groups, '_vs_', ref_group));

pTable_ensembleg = pMatrix_OMM_SPF_vsGF(idxFC,:);
pTable_ensembleg = array2table(pTable_ensembleg, 'RowNames', proteins_ensembleg(idxG,3),...
    'VariableNames', strcat('log2FC_',contrast_groups, '_vs_', ref_group));

fdrTable_ensembleg = fdrMatrix_OMM_SPF_vsGF(idxFC,:);
fdrTable_ensembleg = array2table(fdrTable_ensembleg, 'RowNames', proteins_ensembleg(idxG,3),...
    'VariableNames', strcat('log2FC_',contrast_groups, '_vs_', ref_group));
% save table to file
writetable(fcTable_ensembleg, '.\Output\fcTable_proteins_ensembleg.csv', 'WriteRownames',1);
writetable(pTable_ensembleg, '.\Output\pTable_proteins_ensembleg.csv', 'WriteRownames',1);
writetable(fdrTable_ensembleg, '.\Output\fdrTable_proteins_ensembleg.csv', 'WriteRownames',1);

    