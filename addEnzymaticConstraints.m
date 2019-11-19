function model = addEnzymaticConstraints(model,rxnID,kcat,molMass)
% Add an enzymatic constraint for an existent stoichiometric reaction in
% the form of e_i = v_i*kcat_i
%
% INPUTS
%   model:      Stoichiometric or protein allocation model
%   rxnID:      Character array or cell array containing reaction IDs of
%               reactions for which enzymatic constraints are to be added.
%   kcat:       vector of turnover numbers or kcat values (unit: 1/s)
%   molMasses:  vector conaining the molar masses of enzymes catalyzing the
%               corresponding reactions (unit: g/mol).
% 
% OUTPUTS
%   model:      model with added enzymatic constraints
%
% Author: Tobias B. Alter
% NOV 15, 2019


%% check if protein allocation option exist
if isfield(model,'opt_pa')
    % load options
    opt     = model.opt_pa;
    % check for existent active enzyme constraints
    if ~isfield(opt,'activeEnzymesSector')
        opt.activeEnzymesSector  = []; 
    end
else
    opt     = [];
    opt.activeEnzymesSector  = [];
end

%% check for an total protein constraint
if isfield(opt,'totalProtein')
    if isfield(opt.totalProtein,'cnstrID')
        totalProteinConstraint  = opt.totalProtein.cnstrID;
    else
        totalProteinConstraint  = [];
    end
else
    totalProteinConstraint  = [];
end

%% add enzymatic constraints
if ischar(rxnID)
    rxnID   = {rxnID};
end
% preallocate variables
cnstrIDs        = {};
rxns2ECrxns     = {};
kcat2ECrxns     = [];
ECrxns          = {};
molMass2ECrnxs  = [];
for i=1:length(rxnID)
    % check identifier and kcat value
    rxnPos  = find(strcmp(model.rxns,rxnID{i}));
    if isempty(rxnPos)
        warning(['Reaction "',rxnID{i},'" not found in the model. Skip for active enzyme sector']);
        continue;
    end
    if kcat(i)<0
        warning(['Turnover number for reaction "',rxnID{i},'" is invalid. Skip for active enzyme sector']);
        continue;
    end
    % add active enzyme constraits (as a metabolite)
    constraintID        = ['EA_',rxnID{i}];    % ID in mets vector
    cnstrIDs{end+1,1}   = constraintID;
    model  = addMetabolite(model,constraintID,...
                'metName',['EA_',model.rxnNames{rxnPos}],...
                'b',0,'csense','E'); 
            
    % link metabolic flux and enzyme concentration
    model.S(find(strcmp(model.mets,constraintID)),rxnPos)    = 1;
    enzID     = ['EAR_',model.rxns{rxnPos}];  % enzye ID in reaction vector                           
    model  = addReaction(model,enzID,...
            'reactionName',['EAR_',model.rxnNames{rxnPos}],...
            'metaboliteList',{constraintID},...
            'stoichCoeffList',-(kcat(i)*3600)*1e-6,...  % unit of enzyme concentration is nmol/g
            'reversible',0,...
            'lowerBound',0,...
            'upperBound',1e6,...
            'objectiveCoef',0,...
            'subSystem','Active Enzyme Concentration',...
            'printLevel',0);  

    rxns2ECrxns{end+1,1}    = model.rxns{rxnPos};        
    kcat2ECrxns(end+1,1)    = kcat(i);
    ECrxns{end+1,1}         = enzID;  
    
    % add to total protein constraint
    if ~isempty(totalProteinConstraint)
        model.S(strcmp(model.mets,totalProteinConstraint),...
                        strcmp(model.rxns,enzID))   = molMass(i)*1e-06;
    end
    molMass2ECrnxs(end+1,1)     = molMass(i);        
end

% save and/or append new information
if isempty(opt.activeEnzymesSector)
    % create new info storage 
    opt.activeEnzymesSector.rxns2ECrxns     = rxns2ECrxns;
    opt.activeEnzymesSector.kcat2ECrxns     = kcat2ECrxns;
    opt.activeEnzymesSector.ECrxns          = ECrxns;
    opt.activeEnzymesSector.molMass2ECrnxs  = molMass2ECrnxs;
    opt.activeEnzymesSector.cnstrID2ECrxns  = cnstrIDs;
else
    % append information
    opt.activeEnzymesSector.rxns2ECrxns ...
        = [opt.activeEnzymesSector.rxns2ECrxns;rxns2ECrxns];
    opt.activeEnzymesSector.kcat2ECrxns ...
        = [opt.activeEnzymesSector.kcat2ECrxns;kcat2ECrxns];
    opt.activeEnzymesSector.ECrxns ...
        = [opt.activeEnzymesSector.ECrxns;ECrxns];
    opt.activeEnzymesSector.molMass2ECrnxs ...
        = [opt.activeEnzymesSector.molMass2ECrnxs;molMass2ECrnxs];
    opt.activeEnzymesSector.cnstrID2ECrxns ...
        = [opt.activeEnzymesSector.cnstrID2ECrxns;cnstrIDs];
end
    
    
model.opt_pa  = opt;

end