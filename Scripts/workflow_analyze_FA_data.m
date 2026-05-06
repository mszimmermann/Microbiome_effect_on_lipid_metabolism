dataFile = '.\Data\lipidome_data\Combined_Data_MZimmermann_TotalFA_121020_gutcont_bile.xlsx';
figFolder = '.\Figures\';
outputFolder = '.\Output\';

sheetNames = {'March 2020 Data', 'Oct 2020 Bile fluid', 'Oct 2020 gut contents'};

datTable1 = readtable(dataFile, 'Sheet', sheetNames{1});
datTable2 = readtable(dataFile, 'Sheet', sheetNames{2});
datTable3 = readtable(dataFile, 'Sheet', sheetNames{3});

% multiple gut contents by contents
% for i=9:size(datTable3,2)
%     datTable3{:,i} = datTable3{:,i}.*datTable3{:,7};
% end
%     
mergeColumns = union(datTable1.Properties.VariableNames,...
                     datTable2.Properties.VariableNames);%,...
mergeColumns = union(mergeColumns, datTable3.Properties.VariableNames);
                 
t1colmissing = setdiff(mergeColumns, datTable1.Properties.VariableNames);
t2colmissing = setdiff(mergeColumns, datTable2.Properties.VariableNames);
t3colmissing = setdiff(mergeColumns, datTable3.Properties.VariableNames);

datTable1 = [datTable1 array2table(nan(height(datTable1), numel(t1colmissing)), 'VariableNames', t1colmissing)];
datTable2 = [datTable2 array2table(nan(height(datTable2), numel(t2colmissing)), 'VariableNames', t2colmissing)];
datTable3 = [datTable3 array2table(nan(height(datTable3), numel(t3colmissing)), 'VariableNames', t3colmissing)];
datTable = [datTable1; datTable2; datTable3];

sampleType_unique = unique(datTable.Housing);
sampleTissue_unique = unique(datTable.Matrix);
sampleTime_unique = unique(datTable.Time_h_);
% resort tissues according to GI tract
%sampleTissue_unique = sampleTissue_unique([8 15 13 14 12 2 4 3 1 7 5 9 17 16 6]);
sampleTissue_unique = sampleTissue_unique([13 15 14 12 2 4 3 1 7 5 8 9 17 16 6]);

FAcolumns = datTable.Properties.VariableNames(cellfun(@(x)...
                (contains(x, 'FA') & ~contains(x, 'Total')), datTable.Properties.VariableNames));

FAcolumns_mean = nanmean(datTable{:, FAcolumns});
FAcolumns_mean_top = sort(FAcolumns_mean, 'descend');

% keep top 13 fatty acids
FAcolumns = FAcolumns(FAcolumns_mean>=FAcolumns_mean_top(13));

           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform ranova for all tissues and FA
nrep = 5;
ranovaPtime = ones(length(FAcolumns),length(sampleTissue_unique));
ranovaPtime_type = ones(length(FAcolumns),length(sampleTissue_unique));
aucGF = zeros(length(FAcolumns),length(sampleTissue_unique));
aucSPF = zeros(length(FAcolumns),length(sampleTissue_unique));
aucOMM = zeros(length(FAcolumns),length(sampleTissue_unique));
% store all auc per replicate
aucGFall = cell(length(FAcolumns),length(sampleTissue_unique));
aucSPFall = cell(length(FAcolumns),length(sampleTissue_unique));
aucOMMall = cell(length(FAcolumns),length(sampleTissue_unique));
aucPSPFOMMall = zeros(length(FAcolumns),length(sampleTissue_unique));
aucPSPFGFall = zeros(length(FAcolumns),length(sampleTissue_unique));
aucPOMMGFall = zeros(length(FAcolumns),length(sampleTissue_unique));

for fa_i = 1:length(FAcolumns)
    for tissue_i = 1:length(sampleTissue_unique)
        tissuealltime = nan(nrep*length(sampleType_unique),length(sampleTime_unique));
        tissuemeantime = nan(length(sampleType_unique),length(sampleTime_unique));
        tissuemeantype = cell(nrep *length(sampleType_unique),1);
        idx = 1;
        for type_i = 1:length(sampleType_unique)
            curidx = contains(datTable.Housing, sampleType_unique{type_i}) &...
                     contains(datTable.Matrix, sampleTissue_unique{tissue_i});
            curtime = datTable.Time_h_(curidx);
            curdata = datTable{curidx,...
                ismember(datTable.Properties.VariableNames, FAcolumns{fa_i})};
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
            ranovaPtime(fa_i, tissue_i) = ranovatbl.pValue(1);%LB(1);
            ranovaPtime_type(fa_i, tissue_i) = ranovatbl.pValue(2);%LB(2);
        end
        % calculate auc
        tissuemeantime = sum(tissuemeantime,2);
        aucGF(fa_i, tissue_i) = tissuemeantime(ismember(sampleType_unique, 'GF'));
        aucSPF(fa_i, tissue_i) = tissuemeantime(ismember(sampleType_unique, 'SPF'));
        aucOMM(fa_i, tissue_i) = tissuemeantime(ismember(sampleType_unique, 'OMM12'));
        % calculate auc per replicate
        tissuealltime(isnan(tissuealltime))=0;% do not count nans
        tissuealltime = sum(tissuealltime,2);
        aucGFall{fa_i, tissue_i} = tissuealltime(ismember(tissuemeantype, 'GF'));
        aucSPFall{fa_i, tissue_i} = tissuealltime(ismember(tissuemeantype, 'SPF'));
        aucOMMall{fa_i, tissue_i} = tissuealltime(ismember(tissuemeantype, 'OMM12'));

        [~, aucPSPFOMMall(fa_i, tissue_i)] = ttest2(tissuealltime(ismember(tissuemeantype, 'SPF')),...
                tissuealltime(ismember(tissuemeantype, 'OMM12'))); 
        [~, aucPSPFGFall(fa_i, tissue_i)] = ttest2(tissuealltime(ismember(tissuemeantype, 'SPF')),...
                tissuealltime(ismember(tissuemeantype, 'GF'))); 
        [~, aucPOMMGFall(fa_i, tissue_i)] = ttest2(tissuealltime(ismember(tissuemeantype, 'OMM12')),...
                tissuealltime(ismember(tissuemeantype, 'GF'))); 
        
    end
end
ranovaPtime_typeFDR = reshape(mafdr(ranovaPtime_type(:), 'bhfdr',1),size(ranovaPtime_type));            

aucFDRSPFOMMall = aucPSPFOMMall;
aucFDRSPFGFall = aucPSPFGFall;
aucFDROMMGFall = aucPOMMGFall;

% perform adjustment of p-values per tissue (Benjamini-Hochberg)
for i=1:size(ranovaPtime_typeFDR, 2)
    ranovaPtime_typeFDR(:,i) = mafdr(ranovaPtime_type(:,i), 'bhfdr',1);
    aucFDRSPFOMMall(:,i) = mafdr(aucPSPFOMMall(:,i), 'bhfdr',1);
    aucFDRSPFGFall(:,i) = mafdr(aucPSPFGFall(:,i), 'bhfdr',1);
    aucFDROMMGFall(:,i) = mafdr(aucPOMMGFall(:,i), 'bhfdr',1);
end

% plot pvalues
clusterdata = ranovaPtime_typeFDR;
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
    'ranova_pvalues_tissues_extended_pFDRtissue_0_1.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot fold changes
% clusterdata = (aucSPF+0.000000001) ./ (aucGF+0.000000001);%(aucOMM + 0.000000001);%
clusterdata  = zeros(size(aucSPF,1), size(aucSPF,2)*2);
smallnum = 0.000000001;
clusterdata (:,1:2:end) = log2((aucSPF+smallnum) ./ (aucGF+smallnum));% .*(ranovaPtime_typeFDR<=0.1);
clusterdata (:,2:2:end) = log2((aucOMM+smallnum) ./ (aucGF+smallnum));% .*(ranovaPtime_typeFDT<=0.1);

clusterdata(isnan(clusterdata)) = 0;
clusterdata(isinf(clusterdata)) = 0;

clustercols = reshape([strcat(sampleTissue_unique', ' SPF vs GF') ; ...
                       strcat(sampleTissue_unique', ' OMM12 vs GF')], [],1);

                   
clusterrows = FAcolumns;
clusterrows = cellfun(@(x) strrep(x, '_', '-'), clusterrows, 'unif', 0);


cgo = clustergram(clusterdata,...
            'ColumnLabels', clustercols,...
            'RowLabels', clusterrows,...
            'Cluster', 'column',...
            'symmetric', 1,...
            'colormap', redbluecmap);   
% get the clustergram data as a table in the same order
[~, ~, idx] = intersect(cgo.RowLabels, clusterrows, 'stable');
clusterdata_ordered = clusterdata(idx,:);
clusterdata_ordered = array2table(clusterdata_ordered,...
    'VariableNames', clustercols, 'RowNames', clusterrows(idx));
% resort row labels
clusterdata_ordered = clusterdata_ordered(flipud((1:length(clusterrows))'),:);
writetable(clusterdata_ordered, '.\Output\table_heatmap_data1b.csv',...
    'WriteRowNames',1);


C = findall(gcf,'type','ColorBar');
colorTitleHandle = get(C,'Title');
set(colorTitleHandle ,'String','log2(FC AUC)'); 
orient landscape

% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
%     [figFolder 'Fig1B_clustergram_AUCFC_13FA_tissues_no_sort.pdf'])
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
    [figFolder 'Fig1B_clustergram_AUCFC_13FA_tissues_no_sort_FDR_0_1.pdf'])

% alternative visualization to heatmaps
data_spf = clusterdata(:,1:2:end);
clustercols_spf = clustercols(1:2:end);
data_omm = clusterdata(:,2:2:end);
clustercols_omm = clustercols(2:2:end);

%scatter(data_omm(:), data_cvr(:))
select_fa = [2 3 4];
clustercols_omm_short = cellfun(@(x) x(1:strfind(x, 'OMM')), clustercols_omm, 'unif', 0);
clustercols_spf_short = cellfun(@(x) x(1:strfind(x, 'SPF')), clustercols_spf, 'unif', 0);


fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1)
h = heatmap(clustercols_omm_short, clusterrows(select_fa), data_omm(select_fa,:),...
    'colormap', redbluecmap);   
h.CellLabelFormat = '%.2f';
clim([-4 4])
title('OMM vs GF')
subplot(2,1,2)
h = heatmap(clustercols_spf_short, clusterrows(select_fa), data_spf(select_fa,:),...
    'colormap', redbluecmap);  
h.CellLabelFormat = '%.2f';
clim([-4 4])
title('SPF vs GF')
orient landscape
print(fig, '-vector', '-dpdf', '-r600', '-bestfit', ...
    [figFolder 'Fig1B_clustergram_AUCFC_13FA_tissues_no_sort_FDR_0_1_separate_groups_formatted.pdf'])

% print figure data to file
figuredataOMM = array2table(data_omm(select_fa,:),...
                            'VariableNames', clustercols_omm_short,...
                            'RowNames', clusterrows(select_fa));
writetable(figuredataOMM, '.\Output\table_heatmap_data_Figure1b_omm.csv',...
    'WriteRowNames',1);

figuredataSPF = array2table(data_spf(select_fa,:),...
                            'VariableNames', clustercols_spf_short,...
                            'RowNames', clusterrows(select_fa));
writetable(figuredataSPF, '.\Output\table_heatmap_data_Figure1b_spf.csv',...
    'WriteRowNames',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot AUC of specific FA
mycolors = [0 .45 .7;... % blue         --- GF%
           .89  0.62 0;... % orange   --- DC
            .83 .37 0 ...% vermillion ...  --- CV
            ];  


outFileName = 'bar_auc_fdrttest_';  
for fa_i = 1:length(select_fa)
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    set(groot, 'DefaultTextInterpreter', 'none')
    for tissue_i = 1:length(sampleTissue_unique)
        subplot(3,5,tissue_i)
        bardata = [aucGF(select_fa(fa_i),tissue_i),...
            aucOMM(select_fa(fa_i),tissue_i),...
            aucSPF(select_fa(fa_i),tissue_i)];
        % add colors
        hold on
        for i=1:size(bardata,2)
            bar(i, bardata(:,i), 'FaceColor', mycolors(i,:));
        end
     
        curdata = [aucGFall{select_fa(fa_i),tissue_i};...
            aucOMMall{select_fa(fa_i),tissue_i};...
            aucSPFall{select_fa(fa_i),tissue_i}];
        curx = [ones(5,1)+0.2*(rand(5,1)-0.5);...
            2*ones(5,1)+0.2*(rand(5,1)-0.5);...
            3*ones(5,1)+0.2*(rand(5,1)-0.5)];
        plot(curx, curdata, '.k')
        %boxplot(curdata)
        % plot p-values
        yvalue = max(max(curdata));
        if aucFDROMMGFall(select_fa(fa_i),tissue_i)<=0.1
            plot([1,2], [yvalue*1.1,yvalue*1.1], 'k')
            text(1.5,yvalue*1.1,sprintf('%.2f', aucFDROMMGFall(select_fa(fa_i),tissue_i)))
        end
        if aucFDRSPFGFall(select_fa(fa_i),tissue_i)<=0.1
            plot([1,3],[yvalue*1.2,yvalue*1.2], 'k')
            text(2,yvalue*1.2,sprintf('%.2f', aucFDRSPFGFall(select_fa(fa_i),tissue_i)))
        end
        if aucFDRSPFOMMall(select_fa(fa_i),tissue_i)<=0.1
            plot([2,3],[yvalue*1.3,yvalue*1.3], 'k')
            text(2.5,yvalue*1.3,sprintf('%.2f', aucFDRSPFOMMall(select_fa(fa_i),tissue_i)))
        end
        
        title(sampleTissue_unique{tissue_i})
        set(gca, 'xtick', 1:length(sampleType_unique))
        set(gca, 'xticklabels', sampleType_unique)
        xlim([0.5 3.5])
    % save to file
    end
    sgtitle(FAcolumns{select_fa(fa_i)})
    orient landscape
    print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
                [figFolder outFileName FAcolumns{select_fa(fa_i)}]);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot all FA
mycolors = [0 .6 .5;...  % bluish green    --- CV
            0 .45 .7;... % blue            --- DC
            .8 .4 0];  % vermillion      --- GF%plot experimental data on top

outFileName = 'FAprofiles_across_tissues_extended_multiplied_by_dryweight_finalized';
fig = figure('units','normalized','outerposition',[0 0 1 1]);
set(groot, 'DefaultTextInterpreter', 'none')
    
for fa_i = 1:length(FAcolumns)
    for tissue_i = 1:length(sampleTissue_unique)
        for type_i = 1:length(sampleType_unique)
            curidx = contains(datTable.Housing, sampleType_unique{type_i}) &...
                     contains(datTable.Matrix, sampleTissue_unique{tissue_i});
            curtime = datTable.Time_h_(curidx);
            curdata = datTable{curidx,...
                ismember(datTable.Properties.VariableNames, FAcolumns{fa_i})};
            curunit = unique(datTable.Unit(curidx));
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
            ylabel(curunit)
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
    end
    legend(sampleType_unique, 'Location', 'EastOutside')
    % save to file
    orient landscape
    print(gcf, '-painters', '-dpsc', '-r600', '-bestfit', '-append',...
            outFileName);
    clf(fig)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate mean for modelling

faMeanTable = zeros(length(sampleType_unique)*length(sampleTissue_unique)*length(sampleTime_unique),...
                    length(FAcolumns));
faMeanTable_rows = cell(size(faMeanTable,1),1);
for fa_i = 1:length(FAcolumns)
    idx = 1;
    for type_i = 1:length(sampleType_unique)
        for tissue_i = 1:length(sampleTissue_unique)
            curidx = contains(datTable.Housing, sampleType_unique{type_i}) &...
                     contains(datTable.Matrix, sampleTissue_unique{tissue_i});
            curtime = datTable.Time_h_(curidx);
            curdata = datTable{curidx,...
                ismember(datTable.Properties.VariableNames, FAcolumns{fa_i})};
            curunit = unique(datTable.Unit(curidx));
            % calculate mean per time to plot the line
            curmeantime = zeros(size(sampleTime_unique));
            for j=1:length(sampleTime_unique)
                curmeantime(j) = mean(curdata(curtime==sampleTime_unique(j)));
            end
            faMeanTable(idx:idx+length(curmeantime)-1,fa_i) = curmeantime;
            faMeanTable_rows(idx:idx+length(curmeantime)-1) = ...
                arrayfun(@(x) strcat(sampleType_unique{type_i},'_',...
                                     sampleTissue_unique{tissue_i},'_',...
                                     num2str(x)), sampleTime_unique, 'unif', 0);
            idx = idx+length(curmeantime);
        end
    end
end
faMeanTable_table = array2table(faMeanTable, 'VariableNames', FAcolumns);
faMeanTable_table.Condition = faMeanTable_rows;
writetable(faMeanTable_table, 'fa_mean_table_extended_GITcontents_nmol.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select columns corresponding to FA of interest
select_cols = {'MouseIdentifier', 'Time_h_'};
datTable_for_model = datTable(:, select_cols);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a table with fatty acid per tissue as columns, 
% and add time and mouse group as rows
% save data from mice from the same group with the same ID
fa_selected = {'FA16_0Mz275', 'FA16_0Mz301'};
tissues_selected = sampleTissue_unique;
% remove lung, BAT, gallbladder and bile
tissues_selected(cellfun(@(x) contains(lower(x), 'lung'), tissues_selected))=[];
tissues_selected(cellfun(@(x) contains(lower(x), 'bile'), tissues_selected))=[];
tissues_selected(cellfun(@(x) contains(lower(x), 'ibat'), tissues_selected))=[];
tissues_selected(cellfun(@(x) contains(lower(x), 'gall'), tissues_selected))=[];

mouse_id_unique = unique(datTable.MouseIdentifier);

faMouseTable = nan(length(mouse_id_unique),...
    length(fa_selected)*length(tissues_selected));
faMouseTable_time = zeros(length(mouse_id_unique),1);
faMouseTable_group = cell(length(mouse_id_unique),1);

for mouse_i = 1:length(mouse_id_unique)
    curidx = ismember(datTable.MouseIdentifier, mouse_id_unique{mouse_i});
    curtable = datTable(curidx,:);
    [~, table_order, curorder] = intersect(tissues_selected, curtable.Matrix, 'stable');
    for fa_i = 1:length(fa_selected)
        faMouseTable(mouse_i, (fa_i-1)*length(tissues_selected)+table_order) = curtable{curorder, fa_selected{fa_i}};
    end
    faMouseTable_time(mouse_i) = unique(curtable.Time_h_);
    faMouseTable_group{mouse_i} = unique(curtable.Housing);
end
[y,x]=ndgrid(1:length(tissues_selected),1:length(fa_selected));
faMouseTable_columns = strcat(fa_selected(x(:))','_',tissues_selected(y(:)));
% replace gut content with Content_
faMouseTable_columns = cellfun(@(x) strrep(x, 'gut content ', 'Content_'), faMouseTable_columns, 'unif', 0);

faMouseTable_table = array2table(faMouseTable, 'VariableNames', faMouseTable_columns);
faMouseTable_table.Time = faMouseTable_time;
faMouseTable_table.Group = faMouseTable_group;
% sort columns
faMouseTable_table = faMouseTable_table(:, [{'Time'}, {'Group'}, faMouseTable_columns']);
writetable(faMouseTable_table, 'fa_mouse_table_extended_GITcontents_nmol.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create three tables, one each per mouse group, 
% and add time and mouse group as rows
% save data from mice from the same group with the same ID
fa_selected = {'FA16_0Mz275', 'FA16_0Mz301'};
tissues_selected = sampleTissue_unique;
% remove lung, BAT, gallbladder and bile
tissues_selected(cellfun(@(x) contains(lower(x), 'lung'), tissues_selected))=[];
tissues_selected(cellfun(@(x) contains(lower(x), 'bile'), tissues_selected))=[];
tissues_selected(cellfun(@(x) contains(lower(x), 'ibat'), tissues_selected))=[];
tissues_selected(cellfun(@(x) contains(lower(x), 'gall'), tissues_selected))=[];

mouse_id_unique = unique(datTable.MouseIdentifier);

% GF table
faMouseTable_GF = [];
faMouseTable_GF_time = zeros(length(mouse_id_unique),1);
faMouseTable_GF_group = cell(length(mouse_id_unique),1);
% SPF table
faMouseTable_SPF = [];
faMouseTable_SPF_time = zeros(length(mouse_id_unique),1);
faMouseTable_SPF_group = cell(length(mouse_id_unique),1);
% OMM table
faMouseTable_OMM = [];
faMouseTable_OMM_time = zeros(length(mouse_id_unique),1);
faMouseTable_OMM_group = cell(length(mouse_id_unique),1);

mouse_group_unique = unique(datTable.Housing);
for group_i = 1:length(mouse_group_unique)
    mouse_id_unique = unique(datTable.MouseIdentifier(ismember(datTable.Housing, mouse_group_unique{group_i})));
    faMouseTable = nan(length(mouse_id_unique),length(fa_selected)*length(tissues_selected));
    faMouseTable_time = zeros(length(mouse_id_unique),1);
    faMouseTable_group = cell(length(mouse_id_unique),1);
    for mouse_i = 1:length(mouse_id_unique)
        curidx = ismember(datTable.MouseIdentifier, mouse_id_unique{mouse_i});
        curtable = datTable(curidx,:);
        [~, table_order, curorder] = intersect(tissues_selected, curtable.Matrix, 'stable');
        for fa_i = 1:length(fa_selected)
            faMouseTable(mouse_i, (fa_i-1)*length(tissues_selected)+table_order) = curtable{curorder, fa_selected{fa_i}};
        end
        faMouseTable_time(mouse_i) = unique(curtable.Time_h_);
        faMouseTable_group{mouse_i} = unique(curtable.Housing);
    end
    switch mouse_group_unique{group_i}
        case 'GF'
            faMouseTable_GF = faMouseTable;
            faMouseTable_GF_time = faMouseTable_time;
            faMouseTable_GF_group = faMouseTable_group;
        case 'SPF'
            faMouseTable_SPF = faMouseTable;
            faMouseTable_SPF_time = faMouseTable_time;
            faMouseTable_SPF_group = faMouseTable_group;
        case 'OMM12'
            faMouseTable_OMM = faMouseTable;
            faMouseTable_OMM_time = faMouseTable_time;
            faMouseTable_OMM_group = faMouseTable_group;
    end
end
            
            
[y,x]=ndgrid(1:length(tissues_selected),1:length(fa_selected));
faMouseTable_columns = strcat(fa_selected(x(:))','_',tissues_selected(y(:)));
% replace gut content with Content_
faMouseTable_columns = cellfun(@(x) strrep(x, 'gut content ', 'Content_'), faMouseTable_columns, 'unif', 0);

faMouseTable_GF_columns = cellfun(@(x) strrep(x, 'Mz275_', 'Mz275_GF_'), faMouseTable_columns, 'unif', 0);
faMouseTable_GF_columns = cellfun(@(x) strrep(x, 'Mz301_', 'Mz301_GF_'), faMouseTable_GF_columns, 'unif', 0);
faMouseTable_OMM_columns = cellfun(@(x) strrep(x, 'Mz275_', 'Mz275_OMM_'), faMouseTable_columns, 'unif', 0);
faMouseTable_OMM_columns = cellfun(@(x) strrep(x, 'Mz301_', 'Mz301_OMM_'), faMouseTable_OMM_columns, 'unif', 0);
faMouseTable_SPF_columns = cellfun(@(x) strrep(x, 'Mz275_', 'Mz275_SPF_'), faMouseTable_columns, 'unif', 0);
faMouseTable_SPF_columns = cellfun(@(x) strrep(x, 'Mz301_', 'Mz301_SPF_'), faMouseTable_SPF_columns, 'unif', 0);


faMouseTable_GF_table = array2table(faMouseTable_GF, 'VariableNames', faMouseTable_GF_columns);
% separate between mice from the same time point 
curtime = faMouseTable_GF_time*10;
for i=1:length(curtime)
    while nnz(curtime==curtime(i))>1
        curtime(i) = curtime(i)+1;
    end
end
faMouseTable_GF_table.Time = curtime;
%faMouseTable_GF_table.Group = faMouseTable_GF_group;
faMouseTable_GF_table = sortrows(faMouseTable_GF_table, 'Time');

faMouseTable_SPF_table = array2table(faMouseTable_SPF, 'VariableNames', faMouseTable_SPF_columns);
% separate between mice from the same time point 
curtime = faMouseTable_SPF_time*10;
for i=1:length(curtime)
    while nnz(curtime==curtime(i))>1
        curtime(i) = curtime(i)+1;
    end
end
faMouseTable_SPF_table.Time = curtime;
%faMouseTable_SPF_table.Group = faMouseTable_SPF_group;
faMouseTable_SPF_table = sortrows(faMouseTable_SPF_table, 'Time');

faMouseTable_OMM_table = array2table(faMouseTable_OMM, 'VariableNames', faMouseTable_OMM_columns);
% separate between mice from the same time point 
curtime = faMouseTable_OMM_time*10;
for i=1:length(curtime)
    while nnz(curtime==curtime(i))>1
        curtime(i) = curtime(i)+1;
    end
end
faMouseTable_OMM_table.Time = curtime;
%faMouseTable_OMM_table.Group = faMouseTable_OMM_group;
faMouseTable_OMM_table = sortrows(faMouseTable_OMM_table, 'Time');

% join tables for three groups
faMouseTable_table = outerjoin(faMouseTable_GF_table, faMouseTable_OMM_table, 'MergeKeys',true);
faMouseTable_table = outerjoin(faMouseTable_table, faMouseTable_SPF_table, 'MergeKeys',true);
% sort columns
faMouseTable_columns = faMouseTable_table.Properties.VariableNames;
faMouseTable_columns(cellfun(@(x) isequal(x, 'Time'), faMouseTable_columns))=[];
faMouseTable_table = faMouseTable_table(:, [{'Time'}, faMouseTable_columns]);
% convert time back to hours
faMouseTable_table.Time = floor(faMouseTable_table.Time/10);
writetable(faMouseTable_table, 'fa_mouse_table_combined_for_extended_model_nmol.csv');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save ranova results to file
fid = fopen([outputFolder, 'labeled_lipids_ranova_auc_foldchange_results.csv'], 'w');
aucFC_SPF_GF = log2(aucSPF ./ aucGF);
aucFC_OMM_GF = log2(aucOMM ./ aucGF);

% first print AUC values
fprintf(fid, 'Fatty acid');
for i=1:length(sampleTissue_unique)
    fprintf(fid, ',"AUC_SPF_%s","AUC_OMM_%s","AUC_GF_%s","FC_SPF_GF_%s","FC_OMM_GF_%s","ANOVA_PF_%s","ANOVA_FDR_%s"',...
                sampleTissue_unique{i},sampleTissue_unique{i},sampleTissue_unique{i},sampleTissue_unique{i},...
                sampleTissue_unique{i},sampleTissue_unique{i},sampleTissue_unique{i});
end
fprintf(fid, '\n');
for i=1:length(FAcolumns)
    fprintf(fid, '"%s"', FAcolumns{i});
    for j=1:length(sampleTissue_unique)
        fprintf(fid, ',%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f',...
            aucSPF(i,j), aucOMM(i,j), aucGF(i,j),...
            aucFC_SPF_GF(i,j), aucFC_OMM_GF(i,j),...
            ranovaPtime_type(i,j), ranovaPtime_typeFDR(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);


