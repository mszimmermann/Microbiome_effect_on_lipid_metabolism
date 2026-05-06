% compare liver proteomics with public datasets
% read liver proteomics translated to ensembleg
fcTable_ensembleg = readtable('.\Output\fcTable_proteins_ensembleg.csv', 'ReadRowNames',1);
pTable_ensembleg = readtable('.\Output\pTable_proteins_ensembleg.csv', 'ReadRowNames',1);
fdrTable_ensembleg = readtable('.\Output\fdrTable_proteins_ensembleg.csv', 'ReadRowNames',1);

% remove dots from row names
rowNames = cellfun(@(x) x(1:strfind(x, '.')-1), fcTable_ensembleg.Properties.RowNames, 'unif', 0);
% restore names without .
empty_gene_names = cellfun(@(x) isempty(x), rowNames);
rowNames(empty_gene_names) = fcTable_ensembleg.Properties.RowNames(empty_gene_names);

fcTable_ensembleg.Properties.RowNames = rowNames;

% read joint expression atlas datasets
fcTable_expression_atlas = readtable('.\Output\expression_atlas\totalFCtable_expression_atlas.csv', 'ReadRowNames',1);
totalContrasts_table = readtable('.\Output\expression_atlas\totalContrasts_expression_atlas.csv');

% read expression atlas data from affymetrix datasets
fcTable_expression_affy = readtable('.\Output\expression_atlas\totalFCtable_expression_atlas_affymetrix.csv', 'ReadRowNames',1);
totalContrasts_affy_table = readtable('.\Output\expression_atlas\totalContrasts_expression_atlas_affymetrix.csv', 'delim', ',' );


% check if some p-values are recorded as strings and not number
pTable_file = '.\Output\expression_atlas\totalPtable_expression_atlas.csv';
opts = detectImportOptions(pTable_file);
opts.VariableTypes(2:end) = {'double'};
pTable_expression_atlas = readtable(pTable_file, opts, 'ReadRowNames',1);
pTable_expression_atlas.Row = [];

% define p and fc thresholds
pFDRthreshold = 0.05;
fcThreshold = log2(1.5);

fcTable_expression_atlas_pvaluefiltered = fcTable_expression_atlas;
fcTable_expression_atlas_pvaluefiltered{:,:} = fcTable_expression_atlas_pvaluefiltered{:,:}.*(pTable_expression_atlas{:,:}<=pFDRthreshold);


%join tables by protein ids
fcTable_joint = innerjoin(fcTable_ensembleg, fcTable_expression_atlas, 'Keys','Row');
%fcTable_joint = innerjoin(fcTable_ensembleg, fcTable_expression_atlas_pvaluefiltered, 'Keys','Row');
% remove columns that are called p_value (should not be needed, but for
% testing it happened that the columns were mixed
remove_cols = cellfun(@(x) contains(x, 'p_value'), fcTable_joint.Properties.VariableNames);
fcTable_joint(:, remove_cols) = [];

% add affymetrix data
fcTable_joint = innerjoin(fcTable_joint, fcTable_expression_affy, 'Keys','Row');
fcTable_joint_Contrasts = [totalContrasts_table; totalContrasts_affy_table];


% run PCA analysis
fcValues = fcTable_joint{:,:}';
fcValues(:, isnan(sum(fcValues))) = [];

testSTD = std(fcValues, [], 2);
histogram(testSTD);
fcValues(testSTD>2,:) = [];


%fcValues = zscore(fcValues, 0, 2);
%boxplot(fcValues')

[coeff,score,latent,tsquared,explained] = pca(fcValues);
explained(1:5)
figure
scatter(score(:,1), score(:,2))
hold on
scatter(score(1,1), score(1,2),'r')
scatter(score(2,1), score(2,2),'b')

fcValues_dist = squareform(pdist(fcValues, 'cityblock'));
%clustergram(fcValues_dist)

bar(fcValues_dist(:,1))
scatter(fcValues_dist(:,1), fcValues_dist(:,2))

% sort by distance
[sorted_values, sorted_idx] = sort(fcValues_dist(:,1), 'ascend');

test100 = sorted_idx(1:sorted_values>=50);
test100_ids = fcTable_joint.Properties.VariableNames(test100)'; 
test100_ids = cellfun(@(x) strrep(x, '_log2foldchange', ''), test100_ids, 'unif', 0);
test100_ids(cellfun(@(x) contains(x, 'log2FC'), test100_ids))=[];

[~, ~, idx] = intersect(test100_ids, totalContrasts_table.ContrastID, 'stable');
test100_ids(:,2) = table2cell(totalContrasts_table(idx, 'ContrastName'));

% find MYD88
gene_interest_name = 'myd88';%'nfkb';%
myd88_id = cellfun(@(x) contains(lower(x),gene_interest_name), totalContrasts_table.ContrastName);
myd88_contrast_id = table2cell(totalContrasts_table(myd88_id, 'ContrastID'));
myd88_contrast_id = cellfun(@(x) x{1}, myd88_contrast_id, 'unif', 0);
myd88_column_names = cellfun(@(x) strcat(x, '_log2foldchange'), myd88_contrast_id, 'unif', 0);
fcValues_dist(ismember(fcTable_joint.Properties.VariableNames, myd88_column_names),1)

% find overlap of significant ids
% select proteins based on thresholds in both groups

select_proteins = ((abs(fcTable_ensembleg.log2FC_OMM12_vs_GF)>=fcThreshold) &...
                   (fdrTable_ensembleg.log2FC_OMM12_vs_GF<=fcThreshold) &...
                   (abs(fcTable_ensembleg.log2FC_SPF_vs_GF)>=fcThreshold) &...
                   (fdrTable_ensembleg.log2FC_SPF_vs_GF<=fcThreshold));

select_proteins_ids = fcTable_ensembleg.Properties.RowNames(select_proteins);
select_proteins_table = fcTable_joint(select_proteins_ids,:);
select_proteins_table_values = select_proteins_table{:,:};
%boxplot(select_proteins_table_values)
testFCnum = sum(abs(select_proteins_table_values)>=fcThreshold); 

testFCnum(ismember(fcTable_joint.Properties.VariableNames, myd88_column_names))

testFCnumtotal = sum(abs(fcValues)>=fcThreshold, 2); 

testFCnum = sum(abs(select_proteins_table_values)>=fcThreshold)./(testFCnumtotal'); 
testFCnum(ismember(fcTable_joint.Properties.VariableNames, myd88_column_names))

testFCnumtotal(ismember(fcTable_joint.Properties.VariableNames, myd88_column_names))

[sorted_values, sorted_idx] = sort(testFCnum, 'descend');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get data from Duparc
dataDuparc = readtable('.\Data\Duparc_et_al\GSE73489_series_matrix_dataonly.txt');
affyIDmap = readtable('.\Data\Duparc_et_al\GSE73489_probe_ids');
%map probe id to ensemble gene id
[test, idxDuparc, idxMap] = intersect(dataDuparc.ID_REF, affyIDmap.affy_mogene_2_1_st_v1);
dataDuparc_ensembl = dataDuparc(idxDuparc,:);
dataDuparc_ensembl.('ensembl_gene_id') = table2cell(affyIDmap(idxMap, 'ensembl_gene_id'));
dataDuparc_ensembl.('gene_name') = table2cell(affyIDmap(idxMap, 'external_gene_name'));
% calculate fold change between groups:
% "GSM1895976" "Liver MyD88 wild-type mice fed with control diet" (g1)
% "GSM1895977"	"Liver MyD88 KO mice fed with control diet" (g2)
% "GSM1895978"	"Liver MyD88 WT fed with high-fat diet"	(g3)
% "GSM1895979"  "Liver MyD88 KO fed with high-fat diet" (g4)
dataDuparc_ensembl.('GSE73489_g2_g1_log2foldchange') = dataDuparc_ensembl{:,3} - dataDuparc_ensembl{:,2};
dataDuparc_ensembl.('GSE73489_g3_g1_log2foldchange') = dataDuparc_ensembl{:,4} - dataDuparc_ensembl{:,2};
dataDuparc_ensembl.('GSE73489_g4_g1_log2foldchange') = dataDuparc_ensembl{:,5} - dataDuparc_ensembl{:,2};
% leave max fold change per gene
dataDuparc_ensembl_maxFC = dataDuparc_ensembl;
[~, idx] = unique(dataDuparc_ensembl_maxFC.ensembl_gene_id);
dataDuparc_ensembl_maxFC = dataDuparc_ensembl_maxFC(idx,:);
for i=1:length(idx)
    cur_id = dataDuparc_ensembl_maxFC{i, 'ensembl_gene_id'};
    curfc = dataDuparc_ensembl(ismember(dataDuparc_ensembl.ensembl_gene_id, cur_id),:);
    if height(curfc)>1
        [~, maxidx] = max(abs(curfc.GSE73489_g2_g1_log2foldchange));
        dataDuparc_ensembl_maxFC(i,:) = curfc(maxidx,:);        
    end
end
dataDuparc_ensembl_maxFC.Properties.RowNames = dataDuparc_ensembl_maxFC.ensembl_gene_id;

% save descriptions of comparisons
ContrastID = {'GSE73489_g2_g1'; 'GSE73489_g3_g1'; 'GSE73489_g4_g1'};
ContrastName = {'Liver MyD88 KO mice vs WT fed with control diet';...
                'Liver WT fed with high-fat diet vs WT control diet';...
                'Liver MyD88 KO fed with high-fat diet vs WT control diet'};
dataDuparc_description = table(ContrastID, ContrastName);

% get only fold change table and format in the same way as expression atlas
dataDuparc_ensembl_maxFC = dataDuparc_ensembl_maxFC(:,...
    cellfun(@(x) contains(x, '_log2foldchange'), dataDuparc_ensembl_maxFC.Properties.VariableNames));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add Duparc data to the joint table
% leave only those genes that are in the protein table
dataDuparc_ensembl_maxFC = dataDuparc_ensembl_maxFC(...
    intersect(dataDuparc_ensembl_maxFC.Properties.RowNames, fcTable_joint.Properties.RowNames),:);
% join setting missing values in Duparc table to nan
fcTable_joint = outerjoin(fcTable_joint, dataDuparc_ensembl_maxFC, 'Keys', 'Row');
fcTable_joint_Contrasts = [fcTable_joint_Contrasts; dataDuparc_description];
fcTable_joint_Contrasts.ContrastID = cellfun(@(x) strrep(x, '-', '_'), fcTable_joint_Contrasts.ContrastID, 'unif', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate overrepresentation of changing genes in the protein data
pFDRthreshold = 0.05;
fcThreshold = log2(1.5);

select_proteins = ((abs(fcTable_ensembleg.log2FC_OMM12_vs_GF)>=fcThreshold) &...
                   (fdrTable_ensembleg.log2FC_OMM12_vs_GF<=fcThreshold) &...
                   (abs(fcTable_ensembleg.log2FC_SPF_vs_GF)>=fcThreshold) &...
                   (fdrTable_ensembleg.log2FC_SPF_vs_GF<=fcThreshold));
select_proteins_ID = fcTable_ensembleg.Properties.RowNames(select_proteins);
fisherP = zeros(width(fcTable_joint),1);
fisherOddsRatio = zeros(width(fcTable_joint),1);
fisherConfInt = cell(width(fcTable_joint),1);
fisherOverlap = zeros(width(fcTable_joint),1);
fisherGroupSignificant = zeros(width(fcTable_joint),1);

for i=1:width(fcTable_joint)
    % identify proteins changing in the current column
    cur_select_proteins = (abs(fcTable_joint{:,i})>=fcThreshold);
    cur_select_proteins_ID = fcTable_joint.Properties.RowNames(cur_select_proteins);
    cur_intersect_ID = (intersect(cur_select_proteins_ID, select_proteins_ID));
    cur_intersect = length(cur_intersect_ID);
    cur_select_nochange = length(setdiff(select_proteins_ID, cur_select_proteins_ID));
    cur_noselect_change = length(setdiff(cur_select_proteins_ID, select_proteins_ID));
    cur_noselect_nochange = length(cur_select_proteins)-cur_select_nochange-...
                                cur_noselect_change-cur_intersect;
    x = table([[cur_intersect;cur_select_nochange],[cur_noselect_change;cur_noselect_nochange]]);

    [h,p,stats] = fishertest(x);
    fisherP(i) = p;
    fisherOddsRatio(i) = stats.OddsRatio;
    fisherConfInt{i} = stats.ConfidenceInterval;
    fisherOverlap(i) = cur_intersect;
    fisherGroupSignificant(i) = length(cur_select_proteins_ID);
end
% correction for multiple hypothesis testing
fisherFDR = mafdr(fisherP, 'bhfdr', 1);

% sort by FDR
[sorted_values, sorted_idx] = sort(fisherFDR, 'ascend');

test100_ids = fcTable_joint.Properties.VariableNames(sorted_idx)'; 
test100_ids = cellfun(@(x) strrep(x, '_log2foldchange', ''), test100_ids, 'unif', 0);
test100_ids(cellfun(@(x) contains(x, 'log2FC'), test100_ids))=[];

% add description of contrasts
[~, idxtest100, idx] = intersect(test100_ids, fcTable_joint_Contrasts.ContrastID, 'stable');
test100_ids(idxtest100,2) = table2cell(fcTable_joint_Contrasts(idx, 'ContrastName'));
test100_table = cell2table(test100_ids, 'VariableNames', {'Contrast_ID', 'Contrast_description'});
test100_table.("FDR_Fisher_test") = sorted_values(3:end);
test100_table.("Gene_overlap") = fisherOverlap(sorted_idx(3:end));
test100_table.("DEGenes_in_the_study") = fisherGroupSignificant(sorted_idx(3:end));
test100_table.("OddsRatio") = fisherOddsRatio(sorted_idx(3:end));
test100_table.("ConfidenceInterval") = fisherConfInt(sorted_idx(3:end));
% save to file
writetable(test100_table, 'C:\Users\mazimmer\Downloads\JosefPaperRevisions2025\liver_proteome_study_comparison_with_affymetrix.csv');

plot_values = sorted_values(3:end);
plot_values = -log10(plot_values);
plot_values(plot_values<0)=0.001;

figure
plot(plot_values, 'LineWidth', 2)
set(gca, 'YScale', 'log')
hold on
plot([0, length(plot_values)], [-log10(0.001), -log10(0.001)], '--')
xlabel('Public gene expression studies')
ylabel({'Enrichment FDR (-log10)','of changing gene/protein sets'})
xlim([0, length(plot_values)])
textlabels = {'Insig1', 'REVERB?', 'H-Ras', 'Myd88', 'Sik1'};
texty = [30 19 13 9 6];
for i=1:length(textlabels)
    text(i*15, texty(i), textlabels{i});
end
axis square
title('Datasets with similar gene expression profiles')
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
    '.\Figures\liver_proteome_study_comparison.pdf');


% make data table used for plotting

% plot p against odds ratio
plotx = fisherOddsRatio(3:end);
ploty = -log10(fisherFDR(3:end));
ngenes = fisherOverlap(3:end);
figure
plot(plotx, ploty, '.', 'Color', [.5 .5 .5])
hold on
select_values = (ploty>-log10(0.001)) & (plotx>10) & (ngenes>10);
plot(plotx(select_values), ploty(select_values), '.', 'MarkerSize', 15, 'Color', [213 94 0]/255)

plotteddatatable = array2table([plotx, ploty],...
    'VariableNames', {'Fischer_Odds_ratio', 'Fischer_FDR_neglog10'},...
    'RowNames', plotlabels);
plotteddatatable.FullText = plotlabels_fulltext;
plotteddatatable.highlight_values = select_values;

writetable(plotteddatatable, '.\Output\table_source_data_overrepresentation_Figure8a.csv',...
    'WriteRowNames',1);

plotlabels = cellfun(@(x) strrep(x, '_log2foldchange', ''), fcTable_joint.Properties.VariableNames(3:end), 'unif', 0);
plotlabels = plotlabels(select_values);
[~, idxplot, idx] = intersect(plotlabels, fcTable_joint_Contrasts.ContrastID, 'stable');
plotlabels_fulltext = plotlabels';
plotlabels_fulltext(idxplot) = fcTable_joint_Contrasts.ContrastName(idx);
plotlabels_genenames = {'Chromogranin A', 'Glycogen synthase', 'Glucocorticoid receptor', 'Myd88', 'Mdr2', 'CAR', 'Nr1d1', 'IGF-1',...
    'Rassf1a', 'Sav1', 'Fxr', 'SGLT5', 'JAK2', 'Pit1', 'HDAC3'};
plotlabels_genesonly = cell(size(plotlabels_fulltext));
for i=1:length(plotlabels_fulltext)
    found_genes = cellfun(@(x) contains(lower(plotlabels_fulltext{i}), lower(x)), plotlabels_genenames);
    plotlabels_genesonly{i} = strjoin(plotlabels_genenames(found_genes));
end
% add labels to the plot
plotx = plotx(select_values);
ploty = ploty(select_values);
for i=1:length(plotlabels_genesonly)
    text(plotx(i), ploty(i), plotlabels_genesonly{i})
end
axis square
xlabel('Odds ratio')
ylabel('adjusted p-value, -log10')
legend({'Public mouse studies, RNAseq', 'Public studies with the most significant gene set overlap'}, ...
    'Location', 'SouthOutside')
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
    '.\Figures\liver_proteome_study_comparison_volcano_affymetrix.pdf');

plotgenes_unique = unique(plotlabels_genesonly);
plotgenes_unique_count = cellfun(@(x) nnz(ismember(plotlabels_genesonly, x)), plotgenes_unique);
plotgenes_unique = table(plotgenes_unique, plotgenes_unique_count);