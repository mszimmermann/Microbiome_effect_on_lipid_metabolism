%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add subfolders to path
addpath(genpath(['.' filesep]));
figFolder = ['.' filesep 'Figures' filesep];
outputFolder = ['.' filesep 'Output' filesep];

% define files containing the models, data, and tissue weights:
dataFolder = ['.' filesep 'Data' filesep];
modelFolder = ['.' filesep 'Models' filesep];
modelfilenames = 'model_FA_3enterocytes_2FA_absileum_eiWAT_combined_ColonNull.csv';
datafilenames = 'fa_mouse_table_combined_for_extended_model_nmol.csv';
weightfilenames = 'mouse_organ_weights_estimates.csv';

% build the models
    modelfilename = [modelFolder modelfilenames];
    datafilename = [dataFolder datafilenames];
    weightfilename = [dataFolder weightfilenames];
    % define model output files with out_ preffix and input file name
    outfilename = [dataFolder 'out_' datafilenames];
    
    % build model according to the definition in the table 
    [modelGutUniversal] = create_model_from_file(modelfilename);
    % load experimental data
    [t,t_amount, t_mean_amount, metNamesMap, gd, useForFitting, dataVolumes] = ...
        load_data_from_combined_table(datafilename, weightfilename);
    
    
    % scale the values by 1000 
    % (updated: no need to scale anymore since weights are adjusted in the weight file)
    scale_factor = 1;
    t_Time = t_amount.Time;
    t_mean_Time = t_mean_amount.Time;
    
    t_amount.Variables = t_amount.Variables*scale_factor;
    t_mean_amount.Variables = t_mean_amount.Variables*scale_factor;
    t_amount.Time = t_Time;
    t_mean_amount.Time = t_mean_Time;
    
    % scale input to the model as well (from the model file)
    for i=1:length(modelGutUniversal.Species)
        modelGutUniversal.Species(i).Value = ...
            modelGutUniversal.Species(i).Value*scale_factor;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract model parameter names and initial values
    parameterNames = cell(size(modelGutUniversal.Parameters));
    parameterInitialValues = zeros(size(modelGutUniversal.Parameters));
    parameterNames_fitdata = zeros(size(modelGutUniversal.Parameters));
    for i=1:length(modelGutUniversal.Parameters)
        parameterNames{i} = modelGutUniversal.Parameters(i).Name;
        parameterNames_fitdata(i) = str2double(modelGutUniversal.Parameters(i).Notes);
        parameterInitialValues(i) = modelGutUniversal.Parameters(i).Value;
        if parameterNames_fitdata(i)~=0
            modelGutUniversal.Parameters(i).Value = 0;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIMULATE MODEL WITH data from datafile to estimate host parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %INITIAL PARAMETER INITIALIZATION 
    initPar = cell(nnz(parameterNames_fitdata==1),1);
    initParValue = zeros(nnz(parameterNames_fitdata==1),1);
    idx =1;
    for i=1:length(parameterNames)
        if parameterNames_fitdata(i)==1
           initPar{idx} = ['log(' parameterNames{i}, ')'];
           initParValue(idx) = parameterInitialValues(i);
           idx = idx+1;
        end
    end
    initParValue(initParValue==0) = rand(nnz(initParValue==0),1);
      
    % Possibly run the model several times to select the best fit
    nruns = 1; %number of runs to run the model
    randomRuns_init = zeros(nruns, length(initPar));
    randomRuns_results = cell(nruns,1);

    for run_i = 1:nruns
        init = rand(1,length(initPar));
        for i=1:length(init)
            if initParValue(i)
                init(i) = initParValue(i);
            end
        end
        randomRuns_init(run_i,:) = init;
    end
    for run_i = 1:nruns
        init = randomRuns_init(run_i,:);
        estimated_parameters = estimatedInfo(initPar, 'InitialValue',init);%'
        responseMap = strcat(metNamesMap(useForFitting,1), ' = ',...
                         metNamesMap(useForFitting,2));

        % optimization step
        rng('default')
        globalMethod = 'ga';
        options = optimoptions(globalMethod);
        hybridMethod = 'fminsearch';
        hybridopts = optimset('Display', 'none');
        options = optimoptions(options, 'HybridFcn', {hybridMethod, hybridopts});
        
        gd = groupedData(t_mean_amount);
    
        resultsHost = sbiofit(modelGutUniversal, gd,responseMap,estimated_parameters,[],...
                          globalMethod,options,'pooled',true);
                      
        % save fitting results
        randomRuns_results{run_i} = resultsHost;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
    simulate_pbpk_model_combined_separate_CI(t_mean_amount, resultsHost, ...
        metNamesMap, useForFitting)

    suptitle('Model with no group-specific parameters') 
    orient landscape
    
    print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
            [figFolder 'FigSup_general_model_fits.pdf']);

 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get parameters to vary
    parameters_to_vary = cell(size(modelGutUniversal.Parameters));
    for i=1:length(modelGutUniversal.Parameters)
        parameters_to_vary{i} = modelGutUniversal.Parameters(i).Name;
    end
    % leave only the ones starting with k_
    parameters_to_vary(cellfun(@(x) ~isequal(x(1), 'k'), parameters_to_vary))=[];
    models_to_vary = cell(size(parameters_to_vary));
    results_to_vary = cell(size(parameters_to_vary));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for par_i=1:length(parameters_to_vary)
        [modelGutVariedPar] = vary_model_parameter_per_group(modelGutUniversal, ...
                                                parameters_to_vary{par_i});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % extract model parameter names and initial values
        parameterNames = cell(size(modelGutVariedPar.Parameters));
        parameterInitialValues = zeros(size(modelGutVariedPar.Parameters));
        parameterNames_fitdata = zeros(size(modelGutVariedPar.Parameters));
        for i=1:length(modelGutVariedPar.Parameters)
            parameterNames{i} = modelGutVariedPar.Parameters(i).Name;
            parameterNames_fitdata(i) = str2double(modelGutVariedPar.Parameters(i).Notes);
            parameterInitialValues(i) = modelGutVariedPar.Parameters(i).Value;
            if parameterNames_fitdata(i)~=0
                modelGutVariedPar.Parameters(i).Value = 0;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %INITIAL PARAMETER INITIALIZATION 
        initPar = cell(nnz(parameterNames_fitdata==1),1);
        initParValue = zeros(nnz(parameterNames_fitdata==1),1);
        idx =1;
        for i=1:length(parameterNames)
            if parameterNames_fitdata(i)==1
               initPar{idx} = ['log(' parameterNames{i}, ')'];
               initParValue(idx) = parameterInitialValues(i);
               idx = idx+1;
            end
        end
        initParValue(initParValue==0) = rand(nnz(initParValue==0),1);

        % Possibly run the model several times to select the best fit
        init = rand(1,length(initPar));
        for i=1:length(init)
            if initParValue(i)
                init(i) = initParValue(i);
            end
        end
        estimated_parameters = estimatedInfo(initPar, 'InitialValue',init);%'
        responseMap = strcat(metNamesMap(useForFitting,1), ' = ',...
                         metNamesMap(useForFitting,2));

            % optimization step
        rng('default')
        globalMethod = 'ga';
        options = optimoptions(globalMethod);
        hybridMethod = 'fminsearch';
        hybridopts = optimset('Display', 'none');
        options = optimoptions(options, 'HybridFcn', {hybridMethod, hybridopts});
        resultsHost_varied = sbiofit(modelGutVariedPar, gd,responseMap,estimated_parameters,[],...
                          globalMethod,options,'pooled',true);
                          
         
        % save varied models and fitting results
        models_to_vary{par_i} = modelGutVariedPar;
        results_to_vary{par_i} = resultsHost_varied;
        fprintf('Varied the model parameter %s\n', parameters_to_vary{par_i});
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check MSE, AIC and other parameters of results
    results_mse = zeros(length(results_to_vary)+1,1);
    results_sse = zeros(length(results_to_vary)+1,1);
    results_aic = zeros(length(results_to_vary)+1,1);
    results_bic = zeros(length(results_to_vary)+1,1);
    
    results_mse_delta = zeros(length(results_to_vary),1);
    results_aic_delta = zeros(length(results_to_vary),1);
    results_bic_delta = zeros(length(results_to_vary),1);
    
    results_mse(1) = resultsHost.MSE;
    results_sse(1) = resultsHost.SSE;
    results_aic(1) = resultsHost.AIC;
    results_bic(1) = resultsHost.BIC;
      
    for i=1:length(results_to_vary)
        results_mse(i+1) = results_to_vary{i}.MSE;
        results_sse(i+1) = results_to_vary{i}.SSE;
        results_aic(i+1) = results_to_vary{i}.AIC;
        results_bic(i+1) = results_to_vary{i}.BIC;
        
        results_mse_delta(i) = (results_mse(i+1)-results_mse(1))/results_mse(1);
        results_aic_delta(i) = (results_aic(i+1)-results_aic(1));
        results_bic_delta(i) = (results_bic(i+1)-results_bic(1));
    end
    
    figure
    subplot(2,2,1)
    plot(results_mse)
    set(gca, 'XTick', 1:length(parameters_to_vary)+1)
    set(gca, 'XTickLabel', [{'Original'}; parameters_to_vary])
    title('MSE (mean squared error)')

    subplot(2,2,2)
    plot(results_sse)
    set(gca, 'XTick', 1:length(parameters_to_vary)+1)
    set(gca, 'XTickLabel', [{'Original'}; parameters_to_vary])
    title('SSE (sum squared error)')
    
    subplot(2,2,3)
    plot(results_aic)
    set(gca, 'XTick', 1:length(parameters_to_vary)+1)
    set(gca, 'XTickLabel', [{'Original'}; parameters_to_vary])
    title('AIC (Akaike information criterion)')
    
    subplot(2,2,4)
    plot(results_bic)
    set(gca, 'XTick', 1:length(parameters_to_vary)+1)
    set(gca, 'XTickLabel', [{'Original'}; parameters_to_vary])
    title('BIC (Bayesian information criterion)')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    set(groot, 'DefaultTextInterpreter', 'none')
    set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
       
    subplot(1,3,1)
    barh(results_mse_delta*100)
    set(gca, 'YTick', 1:length(parameters_to_vary))
    set(gca, 'YTickLabel', parameters_to_vary)
    title('MSE (mean squared error)')
    xlabel('Relative difference in MSE to the basic model')
    ylim([0.5 length(parameters_to_vary)+0.5])
    axis square

    subplot(1,3,2)
    barh(results_aic_delta)
    set(gca, 'YTick', 1:length(parameters_to_vary))
    set(gca, 'YTickLabel', parameters_to_vary)
    title('AIC (Akaike information criterion)')
    xlabel('Difference in AIC to the basic model') 
    ylim([0.5 length(parameters_to_vary)+0.5])
    axis square
    
    subplot(1,3,3)
    barh(results_bic_delta)
    set(gca, 'YTick', 1:length(parameters_to_vary))
    set(gca, 'YTickLabel', parameters_to_vary)
    title('BIC (Bayesian information criterion)')
    xlabel('Difference in BIC to the basic model') 
    ylim([0.5 length(parameters_to_vary)+0.5])
    axis square
    
    orient landscape
    print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
            [figFolder 'Fig2B_model_comparisons.pdf']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % plot best model
    best_idx = find(results_mse==min(results_mse))-1;
    results_plot = results_to_vary{best_idx};
            
    simulate_pbpk_model_combined_separate_CI(t_mean_amount, results_plot, ...
        metNamesMap, useForFitting)
    suptitle(sprintf('Model with %s specific parameter', parameters_to_vary{best_idx})) 
    orient landscape
    
    print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
            [figFolder 'Fig2C_model_fits.pdf']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot best parameters
    par_est = results_plot.ParameterEstimates.Estimate;
    par_est_std = results_plot.ParameterEstimates.StandardError;
 
    par_est_error = max(par_est_std ./ par_est,[],2);

    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    set(groot, 'DefaultTextInterpreter', 'none')

    subplot(2,1,1)
    b = bar(par_est(par_est_error<1,:));
    hold on
    errorbar(1:nnz(par_est_error<1), par_est(par_est_error<1), par_est_std(par_est_error<1), 'k.')
   
     ylim([10 60])
    set(gca, 'XTick', 1:nnz(par_est_error<1))
    set(gca, 'XTickLabel', results_plot.ParameterEstimates.Name(par_est_error<1))
    
    title('Best model coefficients')
   % axis square
    
    subplot(2,1,2)
    b = bar(par_est(par_est_error<1,:));
    hold on
    errorbar(1:nnz(par_est_error<1), par_est(par_est_error<1), par_est_std(par_est_error<1), 'k.')
   
    ylim([0 5])
    set(gca, 'XTick', 1:nnz(par_est_error<1))
    set(gca, 'XTickLabel', results_plot.ParameterEstimates.Name(par_est_error<1))
    
    % title('Best model coefficients')
    %axis square
    
    orient landscape
    
    print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
        [figFolder 'Fig2D_coef_gf_spf_omm_extended_combined_3enterocytes_2FA_std_lt1_zoomed.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print parameter values to file
% general model
writetable(resultsHost.ParameterEstimates,...
    [outputFolder 'General_model_parameter_estimates.csv']);
% best model with varied parameter
resultsHost_varied = results_to_vary{best_idx};
writetable(resultsHost_varied.ParameterEstimates,...
    [outputFolder sprintf('Model_with_%s_specific_parameter_estimates.csv',...
        parameters_to_vary{best_idx})]);

% print model assessment to file
table_model_assessment = table(['Original'; parameters_to_vary],...
    results_mse, [0; results_mse_delta],...
    results_aic, [0; results_aic_delta],...
    results_bic, [0; results_bic_delta],...
    'VariableNames', {'Varied_parameter', 'MSE', 'MSE_relative_delta_to_original',...
    'AIC', 'AIC_delta_to_original',...
    'BIC', 'BIC_delta_to_original'});
writetable(table_model_assessment,...
    [outputFolder 'table_models_assessment_results.csv']);
