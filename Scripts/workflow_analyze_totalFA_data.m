% add colors for plotting
%addpath('.\Util\distinguishable_colors\')

% first read data from one tissue
fileName = '.\Data\lipidome_data\Total_Lipidome_Combined_Data_Zimmermann_010221.xlsx';
figFolder = '.\Figures\';
outputFolder = '.\Output\';

% Set names of tissues/sheets
tissue_unique = {'Bile', 'Liver', 'eWAT', 'iBAT', 'iWAT', 'Duo', 'Jej', 'Ile', 'Col', 'Plasma'};
time_unique = [1,2,6];

group_comparisons = {'SPF', 'GF';...
                      'OMM12', 'GF'};
for comp_i = 1:size(group_comparisons,1)
    group1 = group_comparisons{comp_i,1};
    group2 = group_comparisons{comp_i,2};

    fcidx = 1;

    fcMatrix_tissue_time = nan(500, length(time_unique)*length(tissue_unique));
    pMatrix_tissue_time = ones(500, length(time_unique)*length(tissue_unique));
    fdrMatrix_tissue_time = ones(500, length(time_unique)*length(tissue_unique));

    jointMatrix_tissue_time = ones(500, 500);
    joint_tissue = cell(1, 500);
    joint_mouse = cell(1, 500);
    joint_time = zeros(1, 500);

    totalMets = [];
    joint_column=1;

    for j=1:length(tissue_unique)
        fileSheet = tissue_unique{j};
        % read table from file
        % read numbers in number format
        %%% updated matlab
        %tissueTable_num = readtable(fileName, 'Sheet', fileSheet);
        % read everything in char format
        opts = detectImportOptions(fileName, 'Sheet', fileSheet);
        opts = setvartype(opts, opts.VariableNames, 'char');  
        opts.DataRange = 'A2'; %always read data startig from second line
        tissueTable = readtable(fileName, opts, 'Sheet', fileSheet);
        % older version
        %tissueTable = readtable(fileName, 'Sheet', fileSheet);
        
        % convert data table to matrix
        dataColumns = tissueTable.Properties.VariableNames;
        dataColumns(cellfun(@(x) (~contains(x,'_')) |...
                                 (contains(x,'IS')), dataColumns))=[];

        tissueTissue = table2cell(tissueTable(1,dataColumns));
        tissueTime = table2cell(tissueTable(3,dataColumns));
        tissueTime = cellfun(@(x) str2double(x), tissueTime);

        tissueMouse = table2cell(tissueTable(2,dataColumns));

        % take data from number format
        tissueData = table2array(tissueTable(6:end, dataColumns));
        tissueData = cellfun(@(x) str2double(x), tissueData);
        %tissueData = table2array(tissueTable_num(6:end, dataColumns));
        % read type names from variable names of the table
        tissueMets = table2cell(tissueTable(6:end,1));

        % quantile normalization
        tissueData_norm = tissueData;
        % replace 0 with nan
        tissueData_norm(tissueData_norm==0)=nan;
        tissueData_norm = quantilenorm(tissueData_norm);
        tissueData_norm(isnan(tissueData_norm)) = 0;



        if joint_column==1
            jointMatrix_tissue_time(1:size(tissueData_norm,1),...
                            1:size(tissueData_norm,2)) = tissueData_norm;
            total_mets = tissueMets;
            joint_tissue(1:size(tissueData_norm,2)) = tissueTissue;
            joint_mouse(1:size(tissueData_norm,2)) = tissueMouse;
            joint_time(1:size(tissueData_norm,2)) = tissueTime;
        else
            new_mets = setdiff(tissueMets, total_mets);
            total_mets = [total_mets; new_mets];
            [~, oldidx,newidx] = intersect(total_mets, tissueMets, 'stable');
            jointMatrix_tissue_time(oldidx,...
                  joint_column:(joint_column+size(tissueData_norm,2)-1)) = ...
                    tissueData_norm(newidx,:);
            joint_tissue(joint_column:(joint_column+size(tissueData_norm,2)-1)) = tissueTissue;
            joint_mouse(joint_column:(joint_column+size(tissueData_norm,2)-1)) = tissueMouse;
            joint_time(joint_column:(joint_column+size(tissueData_norm,2)-1)) = tissueTime;
        end

        %calculate fold change and ttest p-value
        for i=1:length(time_unique)
            curtime = time_unique(i);

            fcMatrix = (nanmean(tissueData_norm(:, ismember(tissueMouse, group1) &...
                                                   (tissueTime==curtime)),2)./...
                        nanmean(tissueData_norm(:, ismember(tissueMouse, group2) &...
                                                   (tissueTime==curtime)),2));

            [~, pMatrix] = ttest2(log2(tissueData_norm(:, ismember(tissueMouse, group1) &...
                                                   (tissueTime==curtime))),...
                                  log2(tissueData_norm(:, ismember(tissueMouse, group2) &...
                                                   (tissueTime==curtime))),...
                                              'Dim', 2, 'VarType', 'equal');

            pFDR = mafdr(pMatrix, 'bhfdr', 1);

            if (joint_column==1)
                fcMatrix_tissue_time(1:length(pFDR),fcidx) = fcMatrix;
                pMatrix_tissue_time(1:length(pFDR),fcidx) = pMatrix;
                fdrMatrix_tissue_time(1:length(pFDR),fcidx) = pFDR;

            else
                fcMatrix_tissue_time(oldidx,fcidx) = fcMatrix(newidx);
                pMatrix_tissue_time(oldidx,fcidx) = pMatrix(newidx);
                fdrMatrix_tissue_time(oldidx,fcidx) = pFDR(newidx);

            end
            fcidx = fcidx+1;
        end
        joint_column = joint_column+size(tissueData_norm,2);
    end
    fcMatrix_tissue_time(length(total_mets)+1:end,:) = [];
    pMatrix_tissue_time(length(total_mets)+1:end,:) = [];
    fdrMatrix_tissue_time(length(total_mets)+1:end,:) = [];

    jointMatrix_tissue_time(:, joint_column:end) = [];
    jointMatrix_tissue_time(length(total_mets)+1:end,:) = [];
    joint_tissue(joint_column:end) = [];
    joint_mouse(joint_column:end) = [];
    joint_time(joint_column:end) = [];

    % calculate lipid class
    totalMets_class = total_mets;
    for i=1:length(totalMets_class)
        if contains(totalMets_class{i},'-')
            totalMets_class{i} = totalMets_class{i}(1:strfind(totalMets_class{i},'-')-1);
        elseif contains(totalMets_class{i},' ')
            totalMets_class{i} = totalMets_class{i}(1:strfind(totalMets_class{i},' ')-1);
        end
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save ttest results to file
    %fid = fopen([outputFolder, sprintf('total_lipids_ttest_per_timepoint_results_2025_M2019b_%s_%s.csv', group1,group2)], 'w');
     fid = fopen([outputFolder, sprintf('total_lipids_ttest_per_timepoint_results_052026_M2023a_%s_%s.csv', group1,group2)], 'w');

    clusterdata = log2(fcMatrix_tissue_time);
    clusterdataP = pMatrix_tissue_time;
    clusterdataPFDR = fdrMatrix_tissue_time;

    plotcols = strcat(reshape(repmat(tissue_unique,length(time_unique),1),[],1),'_t',...
                       arrayfun(@(x) num2str(x),repmat(time_unique',length(tissue_unique),1)));
    fprintf(fid, 'Lipid,Lipid class');
    for i=1:length(plotcols)
        fprintf(fid, ',"%s_FC_%s_vs_%s","%s_P_%s_vs_%s","%s_pFDR_%s_vs_%s"',...
            plotcols{i}, group1, group2,...
            plotcols{i}, group1, group2,...
            plotcols{i}, group1, group2);
    end
    fprintf(fid, '\n');
    for i=1:length(total_mets)
        fprintf(fid, '"%s","%s"', total_mets{i}, totalMets_class{i});
        for j=1:size(clusterdata,2)
            fprintf(fid, ',%.3f,%.3f,%.3f',clusterdata(i,j), clusterdataP(i,j),clusterdataPFDR(i,j));
        end
        fprintf(fid,'\n');
    end
    fclose(fid);    

end









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
plotdata = log2(fcMatrix_tissue_time).*(fdrMatrix_tissue_time<=0.05);
plotdata(isnan(plotdata))=0;
% figure
% imagesc(plotdata)
% colormap(redbluecmap)
% caxis([-2 2])

plotcols = strcat(reshape(repmat(tissue_unique,length(time_unique),1),[],1),'_t',...
                   arrayfun(@(x) num2str(x),repmat(time_unique',length(tissue_unique),1)));
plotrows = total_mets;
clustdist = 'euclidean';

zerorows = sum(plotdata~=0,2)==0;
plotdata(zerorows,:)=[];
plotrows(zerorows)=[];

cgo = clustergram(plotdata,...
            'RowLabels', plotrows,...
            'ColumnLabels', plotcols,...
            'ColumnPdist',clustdist,...
            'RowPdist', clustdist,...
            'Cluster', 'column',...
            'DisplayRange', 2,...
            'colormap', redbluecmap,...
            'ImputeFun', @knnimpute,...
            'OptimalLeafOrder', 0);

orient landscape
print(gcf, '-painters', '-dpng', '-r600',...
    sprintf('clustergram_fa_tissues_%s_vs_%s.png',group1, group2))

orient landscape
print(gcf, '-painters', '-dpng', '-r600',...
    sprintf('clustergram_fa_tissues_%s_vs_%s_subset9.png',group1, group2))


% plot scatter plots

fig = figure('units','normalized','outerposition',[0 0 1 1]);
set(groot, 'DefaultTextInterpreter', 'none')

for i=1:size(fcMatrix_tissue_time,2)
    curfc = log2(fcMatrix_tissue_time(:,i));
    curfdr = fdrMatrix_tissue_time(:,i);
    
    select_up = (curfc>=1) & (curfdr<=0.05);
    select_down = (curfc<=-1) & (curfdr<=0.05);
   
    curfdr = -log10(curfdr);
    
    subplot(5,6,i)
    
    scatter(curfc, curfdr, 'k.')
    hold on
     
    scatter(curfc(select_up), curfdr(select_up), 'r.')
    scatter(curfc(select_down), curfdr(select_down), 'b.')
    
    xlim([-4 4])
    ylim([0 6])
    
    title(plotcols(i))
    
    if mod(i-1, 6)==0
        ylabel('FDR, -log10')
    end
    if i>4*6
        xlabel('log2FC')
    end
end
suptitle(sprintf('%s vs %s', group1, group2));

orient landscape
print(gcf, '-painters', '-dpng', '-r600',...
    sprintf('volcano_fa_tissues_%s_vs_%s.png',group1, group2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform ranova for all tissues and FA
FAcolumns = total_mets;
sampleTissue_unique = unique(joint_tissue);
sampleType_unique = unique(joint_mouse);
sampleTime_unique = unique(joint_time);

nrep = 5;
ranovaPtime = ones(length(FAcolumns),length(sampleTissue_unique));
ranovaPtime_type = ones(length(FAcolumns),length(sampleTissue_unique));
aucGF = zeros(length(FAcolumns),length(sampleTissue_unique));
aucSPF = zeros(length(FAcolumns),length(sampleTissue_unique));
aucOMM = zeros(length(FAcolumns),length(sampleTissue_unique));

for fa_i = 1:length(FAcolumns)
    for tissue_i = 1:length(sampleTissue_unique)
        tissuealltime = nan(nrep*length(sampleType_unique),length(sampleTime_unique));
        tissuemeantime = nan(length(sampleType_unique),length(sampleTime_unique));
        tissuemeantype = cell(nrep *length(sampleType_unique),1);
        idx = 1;
        for type_i = 1:length(sampleType_unique)
            curidx = contains(joint_mouse, sampleType_unique{type_i}) &...
                     contains(joint_tissue, sampleTissue_unique{tissue_i});
            curtime = joint_time(curidx);
            curdata = jointMatrix_tissue_time(fa_i,curidx);
            % save into time matrix
            for j=1:length(sampleTime_unique)
                tissuealltime(idx:idx+nnz(curtime==sampleTime_unique(j))-1,j) = ...
                    curdata(curtime==sampleTime_unique(j));
                 tissuemeantime(type_i,j) = nanmean(curdata(curtime==sampleTime_unique(j)));
            end
            tissuemeantype(idx:idx+nrep-1) = sampleType_unique(type_i);
           
            idx = idx+nrep;
        end
        if nnz(~isnan(tissuealltime))>20
            t = table(tissuemeantype,...
                      tissuealltime(:,1),...
                      tissuealltime(:,2),...
                      tissuealltime(:,3),...
                      tissuealltime(:,4),...
                      'VariableNames',{'Type','t0','t1','t2','t6'});    
            %rm = fitrm(t,'t0-t6 ~ Type','WithinDesign',sampleTime_unique);
            % fit to 1-6 because otherwise too many nans
            rm = fitrm(t,'t1-t6 ~ Type','WithinDesign',sampleTime_unique(2:end));
            ranovatbl = ranova(rm);
            ranovaPtime(fa_i, tissue_i) = ranovatbl.pValueLB(1);
            ranovaPtime_type(fa_i, tissue_i) = ranovatbl.pValueLB(2);
        end
        % calculate auc
        tissuemeantime = sum(tissuemeantime,2);
        aucGF(fa_i, tissue_i) = tissuemeantime(ismember(sampleType_unique, 'GF'));
        aucSPF(fa_i, tissue_i) = tissuemeantime(ismember(sampleType_unique, 'SPF'));
        aucOMM(fa_i, tissue_i) = tissuemeantime(ismember(sampleType_unique, 'OMM12'));
    end
end
ranovaPtime_typeFDR = reshape(mafdr(ranovaPtime_type(:), 'bhfdr',1),size(ranovaPtime_type));            

test = ranovaPtime_type(:);
test = test(~isnan(test));
testFDR = mafdr(test, 'bhfdr',1);            
[fdr,q,priori,R2] = mafdr(test','Method','polynomial','Showplot',true);
ranovaPtime_typeFDR = reshape(mafdr(ranovaPtime_type(:), 'bhfdr',1),size(ranovaPtime_type));            

% plot pvalues
clusterdata = ranovaPtime_type;
clusterdata(clusterdata>0.1) = 1;
clusterdata = -log10(clusterdata);
clusterdata(isnan(clusterdata)) = 0;

clusterrows = FAcolumns;
clusterrows = cellfun(@(x) strrep(x, '_', '-'), clusterrows, 'unif', 0);

clustergram(clusterdata,...
            'ColumnLabels', sampleTissue_unique,...
            'RowLabels', clusterrows,...
            'symmetric', 0,...
            'colormap', parula)           

C = findall(gcf,'type','ColorBar');
colorTitleHandle = get(C,'Title');
set(colorTitleHandle ,'String','-log10(p-value)'); 
orient landscape
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
%     'ranova_pvalues_tissues_p_0_1.pdf')
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
%     'ranova_pvalues_tissues_all_p.pdf')
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
    'ranova_pvalues_tissues_extended_p_0_1.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot fold changes
clusterdata
clusterdata = (aucSPF+0.000000001) ./ (aucGF+0.000000001);
%clusterdata = (aucOMM+0.000000001) ./ (aucGF+0.000000001);
clusterdata = log2(clusterdata);
clusterdata = clusterdata.*(ranovaPtime_type<=0.1);
%clusterdata = clusterdata.*(ranovaPtime_typeFDR<=0.1);

clusterdata(isnan(clusterdata)) = 0;
clusterdata(isinf(clusterdata)) = 0;

totalMets_class_unique = unique(totalMets_class);
clusterrows = totalMets_class;% FAcolumns;
clusterrows = cellfun(@(x) strrep(x, '_', '-'), clusterrows, 'unif', 0);

[clusterrows, sortidx] = sort(clusterrows);
clusterdata = clusterdata(sortidx,:);

clustercols = sampleTissue_unique;
removetissues = cellfun(@(x) contains(x, 'gall bladder'), clustercols);
clusterdata = clusterdata(:, ~removetissues);
clustercols = clustercols(~removetissues);

clustershow = sum(clusterdata~=0,2)>0;
%clusterdata(abs(clusterdata)<log2(1.5)) = 0;


%cm = struct('Labels',{'CE', 'LPE', 'LPC', 'PS'},...
%     'Colors',mycolors{:});%[0 0 0], [0.35 0.7 0.9],[0 0.6 0.5],[0.8 0.4 0]});        
cluster_labels = unique(clusterrows(clustershow));
%mycolors = distinguishable_colors(length(cluster_labels));
mycolors = hsv(length(cluster_labels));

cm = struct('Labels',cluster_labels,...
    'Colors',mycolors);
for i=1:length(cm)
    cm(i).Colors = mycolors(i,:);
end

clustergram(clusterdata(clustershow,:),...
            'ColumnLabels', clustercols,...
            'RowLabels', clusterrows(clustershow),...
            'symmetric', 1,...
            'DisplayRange', 2,...
            'Cluster', 'row',...
            'colormap', redbluecmap,... 
            'RowLabelsColor',cm,...
            'LabelsWithMarkers',true)          

C = findall(gcf,'type','ColorBar');
colorTitleHandle = get(C,'Title');
%set(colorTitleHandle ,'String','log2(AUCSPF/AUCGF)'); 
set(colorTitleHandle ,'String','log2(AUC OMM12/AUC GF)'); 
orient landscape
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
%     'ranova_pvalues_tissues_p_0_1.pdf')
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
%     'log2fold_changes_SPF_GF_all_FC.pdf')
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
%     'log2fold_changes_SPF_GF_FC_1_5.pdf')
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
%     'log2fold_changes_SPF_GF_FC_1_5_extended_tissues.pdf')
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
%     'clustergram_ranova_FDR_0_1_AUClog2fold_changes_SPF_GF_colored.pdf')
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
%     'clustergram_ranova_FDR_0_1_AUClog2fold_changes_OMM12_GF_colored.pdf')
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
%     'clustergram_ranova_FDR_0_1_AUClog2fold_changes_OMM12_GF_colored_sortedbyclass.pdf')
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
    'clustergram_ranova_FDR_0_1_AUClog2fold_changes_SPF_GF_colored_sortedbyclass.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plottissues = [{'bile fluid'} ...
               {'Duodenum'}    {'Jejunum'}  {'Ileum'}    {'Colon'}...
               {'Plasma'} {'Liver'} {'eWAT'} {'iWAT'} {'iBAT'}];
[~, ~, plotidx] = intersect(plottissues, sampleTissue_unique, 'stable');

plotdata = zeros(size(aucSPF,1), length(plotidx)*2);
smallnum = 0.000000001;
plotdata(:,1:2:end) = log2((aucSPF(:,plotidx)+smallnum) ./ (aucGF(:,plotidx)+smallnum)) .*(ranovaPtime_type(:,plotidx)<=0.1);
plotdata(:,2:2:end) = log2((aucOMM(:,plotidx)+smallnum) ./ (aucGF(:,plotidx)+smallnum)) .*(ranovaPtime_type(:,plotidx)<=0.1);

plotdata(isnan(plotdata)) = 0;
plotdata(isinf(plotdata)) = 0;

totalMets_class_unique = unique(totalMets_class);
clusterrows = totalMets_class;% FAcolumns; %
clusterrows = cellfun(@(x) strrep(x, '_', '-'), clusterrows, 'unif', 0);

[clusterrows, sortidx] = sort(clusterrows);
plotdata = plotdata(sortidx,:);
plotdata_metnames = total_mets(sortidx);

clustercols = reshape([strcat(plottissues, ' SPF vs GF') ; ...
                       strcat(plottissues, ' OMM12 vs GF')], [],1);

clustershow = sum(plotdata~=0,2)>0;

plotdata(plotdata>2)=2;
plotdata(plotdata<-2)=-2;

% add row colors
plotrows = clusterrows(clustershow);
plotrowcolors = zeros(size(plotrows,1),3);
% calculate number of lipids in each class)
plotrownumitems = zeros(size(plotrows,1),1);

for i=1:length(plotrows)
    plotrowcolors(i,:) = cm(find(ismember({cm.Labels}, plotrows{i}))).Colors;
    plotrownumitems(i) = nnz(ismember(plotrows, plotrows{i}));
end


myfig = figure('units','normalized','outerposition',[0 0.05 0.9 0.9]);

subplot(1,2,1)
hold on
for i=1:length(plotrows)
    h(i) = scatter(1, i, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', plotrowcolors(i,:));
end
ylim([0, length(plotrows)])
[~, uniquerows] = unique(plotrows, 'stable');
legend(h(flipud(uniquerows)), strcat(plotrows(flipud(uniquerows)), arrayfun(@(x) [' (n=' num2str(x) ')'], plotrownumitems(flipud(uniquerows)),'unif',0)));
title('Lipid class legend')


subplot(1,2,2)
h = heatmap(plotdata(clustershow,:));
colormap redbluecmap
caxis([-2 2])

h.XDisplayLabels = clustercols;
h.YDisplayData = flipud(h.YDisplayData);     
title('Fold change of AUC between mouse groups')

orient landscape

print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
    [figFolder 'Fig3D_heatmap_clustergram_ranova_FDR_0_1_AUClog2fold_changes_sortedbyclass.pdf'])

%%%%%%%%%%
% plot two mouse groups separately
myfig = figure('units','normalized','outerposition',[0 0.05 0.9 0.9]);

subplot(1,3,1)
hold on
for i=1:length(plotrows)
    h(i) = scatter(1, i, 's', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', plotrowcolors(i,:));
end
ylim([0, length(plotrows)])
[~, uniquerows] = unique(plotrows, 'stable');
legend(h(flipud(uniquerows)), strcat(plotrows(flipud(uniquerows)), arrayfun(@(x) [' (n=' num2str(x) ')'], plotrownumitems(flipud(uniquerows)),'unif',0)));
title('Lipid class legend')


subplot(1,3,3)
h = heatmap(plotdata(clustershow,1:2:end));
colormap redbluecmap
caxis([-2 2])

h.XDisplayLabels = clustercols(1:2:end);
h.YDisplayData = flipud(h.YDisplayData);     
title('Fold change of AUC between mouse groups')

subplot(1,3,2)
h = heatmap(plotdata(clustershow,2:2:end));
colormap redbluecmap
caxis([-2 2])

h.XDisplayLabels = clustercols(2:2:end);
h.YDisplayData = flipud(h.YDisplayData);     
title('Fold change of AUC between mouse groups')

orient landscape

print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
    [figFolder 'Fig3D_heatmap_clustergram_ranova_FDR_0_1_AUClog2fold_changes_sortedbyclass_seperatedgroups.pdf'])

% plot data to file
plotteddataOMM = array2table(plotdata(clustershow,2:2:end),...
    'VariableNames', clustercols(2:2:end),...
    'RowNames', plotdata_metnames(clustershow));     
plotteddataOMM.LipidClass = plotrows;

writetable(plotteddataOMM, '.\Output\table_source_data_heatmap_Figure3b_omm.csv',...
    'WriteRowNames',1);

plotteddataSPF = array2table(plotdata(clustershow,1:2:end),...
    'VariableNames', clustercols(1:2:end),...
    'RowNames', plotdata_metnames(clustershow));     
plotteddataSPF.LipidClass = plotrows;

writetable(plotteddataSPF, '.\Output\table_source_data_heatmap_Figure3b_spf.csv',...
    'WriteRowNames',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save ranova results to file
fid = fopen([outputFolder, 'total_lipids_ranova_auc_foldchange_results2025_m2023a_numericinput.csv'], 'w');
aucFC_SPF_GF = log2(aucSPF ./ aucGF);
aucFC_OMM_GF = log2(aucOMM ./ aucGF);

% first print AUC values
fprintf(fid, 'Lipid,Lipid class');
for i=1:length(sampleTissue_unique)
    fprintf(fid, ',"AUC_SPF_%s","AUC_OMM_%s","AUC_GF_%s","FC_SPF_GF_%s","FC_OMM_GF_%s","ANOVA_PF_%s","ANOVA_FDR_%s"',...
                sampleTissue_unique{i},sampleTissue_unique{i},sampleTissue_unique{i},sampleTissue_unique{i},...
                sampleTissue_unique{i},sampleTissue_unique{i},sampleTissue_unique{i});
end
fprintf(fid, '\n');
for i=1:length(total_mets)
    fprintf(fid, '"%s","%s"', total_mets{i}, totalMets_class{i});
    for j=1:length(sampleTissue_unique)
        fprintf(fid, ',%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f',...
            aucSPF(i,j), aucOMM(i,j), aucGF(i,j),...
            aucFC_SPF_GF(i,j), aucFC_OMM_GF(i,j),...
            ranovaPtime_type(i,j), ranovaPtime_typeFDR(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mycolors = [0 .6 .5;...  % bluish green    --- CV
            0 .45 .7;... % blue            --- DC
            .8 .4 0];  % vermillion      --- GF%plot experimental data on top

outFileName = 'FAprofiles_across_tissues_extended_total_lipids';
fig = figure('units','normalized','outerposition',[0 0 1 1]);
set(groot, 'DefaultTextInterpreter', 'none')
    
for fa_i = 1:length(FAcolumns)
    for tissue_i = 1:length(sampleTissue_unique)
        for type_i = 1:length(sampleType_unique)
            curidx = contains(joint_mouse, sampleType_unique{type_i}) &...
                     contains(joint_tissue, sampleTissue_unique{tissue_i});
            curtime = joint_time(curidx);
            curdata = jointMatrix_tissue_time(fa_i,curidx);

            %curunit = unique(datTable.Unit(curidx));
            % calculate mean per time to plot the line
            curmeantime = zeros(size(sampleTime_unique));
            curstdtime = zeros(size(sampleTime_unique));
            for j=1:length(sampleTime_unique)
                curmeantime(j) = mean(curdata(curtime==sampleTime_unique(j)));
                curstdtime(j) = std(curdata(curtime==sampleTime_unique(j)));
            end
            subplot(3,6,tissue_i)
            hold on
            scatter(curtime, curdata, 'MarkerFaceColor',mycolors(type_i,:),...
                'MarkerEdgeColor',mycolors(type_i,:))
%             plot(sampleTime_unique, curmeantime, ...
%                         'LineWidth', 2,...
%                         'Color', mycolors(type_i,:))
            errorbar(sampleTime_unique, curmeantime, curstdtime,...
                        'LineWidth', 2,...
                        'Color', mycolors(type_i,:))
            title(sampleTissue_unique{tissue_i})
            %ylabel(curunit)
            ylim([0, inf])
            xlim([0, 6])
        end
    end
    suptitle(FAcolumns{fa_i})
    % add legend
    subplot(3,6,tissue_i+1)
    hold on
    for type_i = 1:length(sampleType_unique)
       % plot(sampleTime_unique, curmeantime, ...
       %                 'LineWidth', 2,...
       %                 'Color', mycolors(type_i,:))
        errorbar(sampleTime_unique, curmeantime, curstdtime,...
                        'LineWidth', 2,...
                        'Color', mycolors(type_i,:))
        ylim([0, inf])
        xlim([0, 6])
    end
    legend(sampleType_unique, 'Location', 'EastOutside')
    % save to file
    orient landscape
    print(gcf, '-painters', '-dpsc', '-r600', '-bestfit', '-append',...
            outFileName);
    clf(fig)
end
         
