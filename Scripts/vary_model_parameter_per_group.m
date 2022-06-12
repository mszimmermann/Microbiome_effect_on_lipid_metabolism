function [modelGutVariedPar] = vary_model_parameter_per_group(modelGutUniversal, parName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get all reactions in the input model that include parameter parName
% and set different parameters for all of them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% copy model from the original
modelGutVariedPar = copyobj(modelGutUniversal);

parName_gf = [parName '_gf'];
parName_spf = [parName '_spf'];
parName_omm = [parName '_omm'];

% find units and value of the parameter of interest
parValue = [];
parUnits = '';
model_params = get (modelGutUniversal, 'Parameters');
for i=1:length(model_params)
    if isequal(model_params(i).Name, parName)
        parValue = model_params(i).Value;
        parUnits = model_params(i).Units;
        % set notes to 0 to not fit
        modelGutVariedPar.Parameters(i).Notes = '0';
        break;
    end
end
% if parValue is empty exit with warning that parameter not found
if i==length(model_params)
    if isempty(parValue)
        fprintf('Exited without creating model: parameter %s not found\n',...
            parName);
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% add new parameters to the model specific for each mouse type
addparameter(modelGutVariedPar, parName_gf, parValue, ...
                 'ValueUnits', parUnits,...
                 'Notes', '1');
addparameter(modelGutVariedPar, parName_omm, parValue, ...
                 'ValueUnits', parUnits,...
                 'Notes', '1');
addparameter(modelGutVariedPar, parName_spf, parValue, ...
                 'ValueUnits', parUnits,...
                 'Notes', '1');
             
% get indeces all reactions that include parameter parName
rxn_indeces = zeros(length(modelGutVariedPar.Reactions),1);
for i=1:length(modelGutVariedPar.Reactions)
    if contains(modelGutVariedPar.Reactions(i).ReactionRate, parName)
        %set new kinetic law based on metabolite names
        if contains(lower(modelGutVariedPar.Reactions(i).ReactionRate), 'gf')
            %add GF parameter
            new_rate = strrep(modelGutVariedPar.Reactions(i).ReactionRate,...
                parName, parName_gf);
        elseif contains(lower(modelGutVariedPar.Reactions(i).ReactionRate), 'spf')
            %add SPF parameter
            new_rate = strrep(modelGutVariedPar.Reactions(i).ReactionRate,...
                parName, parName_spf);
        elseif contains(lower(modelGutVariedPar.Reactions(i).ReactionRate), 'omm')
            %add OMM parameter
            new_rate = strrep(modelGutVariedPar.Reactions(i).ReactionRate,...
                parName, parName_omm);
        end
        % set new rate
        set(modelGutVariedPar.Reactions(i), 'ReactionRate', new_rate);
        % describe changes
        fprintf('Old: %s \t %s\nNew: %s \t %s \n\n',...
                modelGutUniversal.Reactions(i).Reaction,...
                modelGutUniversal.Reactions(i).ReactionRate,...
                modelGutVariedPar.Reactions(i).Reaction,...
                new_rate);
        rxn_indeces(i) = 1;
    end
end
        
