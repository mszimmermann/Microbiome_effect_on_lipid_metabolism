%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function simulate_pbpk_model_combined_separate_CI(t, modelresults, ...
    metNamesMap, useForFitting)
%% Simulate
    simData = sbiopredictionci(modelresults, 'Alpha', 0.5);
    % get results
    simData = simData.Results;
     
    
    mycolors_gf = [0 150 255]/255; %blue [0 .6 .5];  % bluish green
    mycolors_spf = [196 89 17]/255; % brown [.8 .4 0];  % vermillion
    mycolors_omm = [255 147 0]/255; %yellowish [0 .45 .7]; % blue

    simData_response = string(simData.Response);
    simData_response = cellfun(@(x) strrep(x, 'unnamed.',''), simData_response, 'unif', 0); 
    simData_response_group = simData_response;
    simData_response = cellfun(@(x) strrep(x, '_GF', ''), simData_response, 'unif', 0);
    simData_response = cellfun(@(x) strrep(x, '_SPF', ''), simData_response, 'unif', 0);
    simData_response = cellfun(@(x) strrep(x, '_OMM', ''), simData_response, 'unif', 0);
    
    simData_response_unique = unique(simData_response);
    simData_response_group_unique = unique(simData_response_group);
    simData_response_wogroup_unique = simData_response_group_unique;
    simData_response_wogroup_unique = cellfun(@(x) strrep(x, '_GF', ''), simData_response_wogroup_unique, 'unif', 0);
    simData_response_wogroup_unique = cellfun(@(x) strrep(x, '_SPF', ''), simData_response_wogroup_unique, 'unif', 0);
    simData_response_wogroup_unique = cellfun(@(x) strrep(x, '_OMM', ''), simData_response_wogroup_unique, 'unif', 0);
    
    simData_group = simData.Group;
    simData_group_unique = unique(simData_group);
    mycolors_group = {mycolors_gf, mycolors_omm, mycolors_spf};
    
    simData_time = simData.Time;
    simData_estimate = simData.Estimate;
    simData_CI = simData.ConfidenceInterval;
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    lineWidth = 0.5;
    for i=1:length(simData_response_group_unique)
        spidx = find(ismember(simData_response_unique, simData_response_wogroup_unique{i}));
        subplot(3,ceil(length(simData_response_unique)/3),spidx)
        
        % plot different groups
        plotidx = (ismember(simData_response_group, simData_response_group_unique{i}));
        if contains(simData_response_group_unique{i}, 'GF')
            mycolors = mycolors_group{1};
        elseif contains(simData_response_group_unique{i}, 'OMM')
            mycolors = mycolors_group{2};
        else
            mycolors = mycolors_group{3};
        end

        hold on
        if ismember(simData_response_group_unique{i}, metNamesMap(useForFitting,1)) % used for fitting
            h1 = plot(simData_time(plotidx),...
                simData_estimate(plotidx),...
                'LineWidth', 2,'Color', mycolors);
        else
            h2 = plot(simData_time(plotidx),...
                simData_estimate(plotidx), '--', 'LineWidth', 2,'Color', mycolors);
        end

        % plot confidence interval
        X=[simData_time(plotidx)',fliplr(simData_time(plotidx)')]; %#create continuous x value array for plotting
        Y=[simData_CI(plotidx,1)',fliplr(simData_CI(plotidx,2)')];              %#create y values for out and then back
        fill(X,Y,mycolors, 'FaceAlpha', 0.5);  
         
        title(simData_response_unique{spidx})
        xlim([0 max(simData_time)])
        axis square
    end
    
    %plot experimental data on top
   for i=1:length(simData_response_group_unique)
        spidx = find(ismember(simData_response_unique, simData_response_wogroup_unique{i}));
        subplot(3,ceil(length(simData_response_unique)/3),spidx)
   
        hold on
        % plot different groups
        % plot different groups
        if contains(simData_response_group_unique{i}, 'GF')
            mycolors = mycolors_group{1};
        elseif contains(simData_response_group_unique{i}, 'OMM')
            mycolors = mycolors_group{2};
        else
            mycolors = mycolors_group{3};
        end

        if nnz(ismember(metNamesMap(:,1), simData_response_group_unique{i}))
            dataMet = metNamesMap(ismember(metNamesMap(:,1), simData_response_group_unique{i}),2);
            scatter(t.Time, t{:, dataMet},...
                    'MarkerFaceColor', mycolors);
        end
 
    end

    
    set(groot, 'DefaultTextInterpreter', 'none')

    
    %legend([h1 h2 h3], {'Model fit', 'Prediction', 'Experimental data'})
    %legend([h1 h3], {'Model fit', 'Experimental data'})
    