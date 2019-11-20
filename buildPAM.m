function model_pa = buildPAM(model,totProtConc,varargin)
% Add protein allocation constraints to account for the total
% condition-dependent proteome in stoichiometric reconstructions
%
% INPUTS:
%     model:        Irrversivle stoichiometric reconstruction
% 
%     totProtConc   Total condition-dependent protein concentration (unit:
%                   g/g_CDW) (0.258 g/g_CDW)
% 
%     varargin      Optional Inputs provided as 'ParameterName', Value pairs.
%                   It is recommended to account for the active enzyme,
%                   translational, and excess enzyme sector to cover the whole
%                   condition-dependent proteome in the protein allocation
%                   model.
% 
%             *actEnzSector  Struct containing data to build the active Enzyme
%                            protein sector (AEPS). Each specified reaction will be
%                            coupled to an enzyme concentration via a turnover number.
%               .rxnID:         Reaction IDs of enzymatic reactions
% 
%               .kcat:          Turnover numbers (or kcat values) for each enzymatic
%                               reactions (unit: 1/s). Must be the same size as rxnID.
% 
%               .molMass:       Molar mass for each considered enzyme (unit: g/mol). Must
%                               be the same size as rxnID
% 
%             *transSector:  Struct containing data to build the translational/ribosomal
%                            protein sector (TPS) (linear dependent on the growth rate)
%               .bmRxnID        Model Identifier of the biomass formation reaction as
%                               a character vector
%               .TPS_0          Amount of protein allocated to the translational sector at
%                               zero growth (unit: g/g_CDW)
%               .TPS_mu         Amount of protein allocated to the translational sector per
%                               unit increase in the growth rate (unit: g h/g_CDW)
%               .molMass        (Optional) Molar mass of a fictional
%                               ribosome. If assigned, computed enzyme
%                               concentration resamble ribosome
%                               concentrations
% 
%             *excEnzSector: Struct containing data to build the excess enzyme
%                            protein sector (EEPS) (linear dependent on the substrate uptake rate)
%               .subsRxnID      Model Identifier of the substrate uptake reaction as
%                               a character vector
%               .subsUptakeMax  Maximal substrate uptake rate (unit: mmol/g_CDW/h)
%               .EEPS_0         Amount of protein allocated to the excess enzyme sector at
%                               zero substrate uptake (unit: g/g_CDW)
%               .molMass        (Optional) Molar mass of a concentration unit of 
%                               the excess enzymes sector.
% 
%             *customSectors Cell array containing information (in form of structs)to build additional,
%                            custom protein allocation sectors that are linearly dependent on a
%                            user-defined model variable
%               {}.name         % name of the protein sector
%               {}.linRxnID     Model identifier of the reaction from which the custom
%                               protein sector is linearly dependent
%               {}.CPS_0        Intercept of the linear function describing protein
%                               allocation of the custom protein sector
%               {}.CPS_s        Slope of the linear function describing protein
%                               allocation of the custom protein sector. Ca have
%                               negative and positive values.
%               {}.molMass      (Optional) Molar mass of a concentration unit of 
%                               the custom sector
%
% OUTPUT
%
%     model_pa:     Protein allocation model
%
% Author: Tobias B. Alter
% NOV 15, 2019


%% Check input 

% active enzyme sector
if any(strcmp(varargin,'actEnzSector'))
    actEnzSector    = varargin{find(strcmp(varargin,'actEnzSector'))+1};
    if ~(isfield(actEnzSector,'rxnID') && isfield(actEnzSector,'kcat') ...
           && isfield(actEnzSector,'molMass'))
       warning('Data to build an ACTIVE ENZYME PROTEIN SECTOR is incomplete! Continued neglecting this sector')
       actEnzSector     = [];
    else
        if ~(length(actEnzSector.rxnID)==length(actEnzSector.kcat) ...
                && length(actEnzSector.rxnID)==length(actEnzSector.molMass))
            warning('Sizes of input vectors for enzyme molar masses and turnover numbers are inconsistent. Continued neglecting the active enzymes protein sector')
            actEnzSector = [];
        end
    end
else
    actEnzSector = [];
end

% translational protein sector
if any(strcmp(varargin,'transSector'))
    transSector    = varargin{find(strcmp(varargin,'transSector'))+1};
    if ~(isfield(transSector,'bmRxnID') && isfield(transSector,'TPS_0') ...
           && isfield(transSector,'TPS_mu'))
       warning('Data to build an TRANSLATIONAL PROTEIN SECTOR is incomplete! Continued neglecting this sector')
       transSector     = [];
    else
        % check provided biomass reaction identifier
        transSector.bmRxnNum    = find(strcmp(model.rxns,transSector.bmRxnID));
        if isempty(transSector.bmRxnNum)
            warning('Provided biomass equation name could not be found in the model. Continued neglecting the translational protein sector.')
            transSector = [];
        end
    end
else
    transSector = [];
end


% excess enzymes protein sector
if any(strcmp(varargin,'excEnzSector'))
    excEnzSector    = varargin{find(strcmp(varargin,'excEnzSector'))+1};
    if ~(isfield(excEnzSector,'subsRxnID') && isfield(excEnzSector,'subsUptakeMax') ...
           && isfield(excEnzSector,'EEPS_0'))
       warning('Data to build an EXCESS ENZYMES PROTEIN SECTOR is incomplete! Continued neglecting this sector')
       excEnzSector     = [];
    else
        % check provided biomass reaction identifier
        excEnzSector.subsRxnNum     = find(strcmp(model.rxns,excEnzSector.subsRxnID));
        if isempty(excEnzSector.subsRxnNum)
            warning('Provided substrate uptake reaction could not be found in the model. Continued neglecting the excess enzymes protein sector.')
            excEnzSector = [];
        end
    end
else
    excEnzSector = [];
end


% custom protein sector
if any(strcmp(varargin,'customSectors'))
    customSectors    = varargin{find(strcmp(varargin,'customSectors'))+1};
    validSectors   = zeros(length(customSectors),1);
    for i=length(customSectors)
        if ~(isfield(customSectors{i},'linRxnID') && isfield(customSectors{i},'CPS_0') ...
           && isfield(customSectors{i},'CPS_s'))
           warning(['Data to build an CUSTOM PROTEIN SECTOR ',num2str(i),' is incomplete! Continued neglecting this sector'])
        else
            % check provided reaction identifier
            customSectors{i}.linRxnNum  = find(strcmp(model.rxns,customSectors{i}.linRxnID));
            if isempty(customSectors{i}.linRxnNum)
                warning(['Provided reaction could not be found in the model. Continued neglecting customprotein sector ',num2str(i),'.'])
            else
                validSectors(i)    = 1;
            end
        end  
    end
    customSectors   = customSectors(validSectors);
else
    customSectors   = [];
end

% check if any sectors were provided
if isempty(actEnzSector) && isempty(transSector) && isempty(excEnzSector) ...
        && isempty(customSectors)
    error('No valid protein sectors provided.')   
end

% model irreversibilty
if any(model.lb<0)
    error('Model contains reversible reactions. Irreversibilty of model reactions has to be ensured.')
end

%% setup model
model_pa                = model;
opt_pa                  = [];
opt_pa.model            = model;    % save stichiometric model

% correct csense of model
if size(model_pa.csense,2)>1
    % transpose csense vector
    model_pa.csense     = model_pa.csense';
end 


%% Add total protein constraint
fprintf('Add total condition-dependent protein constraint ...\n')

constraintID   = 'TPC';     % constraint ID (goes as metabolite ID)

model_pa      = addMetabolite(model_pa,constraintID,...
                'metName','totalProeinConstraint',...
                'b',0,'csense','E'); 
TPC_pos     = find(strcmp(model_pa.mets,constraintID));

% set constant total protein concentration
model_pa.b(TPC_pos)     = totProtConc*1000; % unit mg/g
% save info
opt_pa.posConstraints.TPC           = TPC_pos;
opt_pa.totalProtein.totalProteinConcentration    = totProtConc;
opt_pa.totalProtein.cnstrID    = constraintID;

model_pa.opt_pa     = opt_pa;

%% Add active enzymes protein sector
% c.f. the GECKO formulation of Sanchez et al. 2017
if ~isempty(actEnzSector)
    fprintf('Add active enzymes sector ...\n')
    
    model_pa    = addEnzymaticConstraints(model_pa,actEnzSector.rxnID,...
                    actEnzSector.kcat,actEnzSector.molMass);
    
    % reload protein allocation opotions
    opt_pa  = model_pa.opt_pa;
       
end




%% Add translational protein sector
if ~isempty(transSector)
    fprintf('Add translational sector ...\n')
    % define molar mass of a concentration unit of the translational sector
    transMolMass_default    = 4.0590394e05; % default E. coli ribosome molar mass [g/mol]
    if ~isfield(transSector,'molMass')
        transMolMass    = transMolMass_default;
    else
        if ~isempty(transSector.molMass)
            transMolMass    = transSector.molMass;
        else
            transMolMass    = transMolMass_default;
        end
    end
    % add constraint
    constraintID    = 'TPS';    % ID in mets vector
    model_pa  = addMetabolite(model_pa,constraintID,...
                'metName','translationalSector',...
                'b',transSector.TPS_0*1000,'csense','E');                  
    transSectConstraint_pos     = find(strcmp(model_pa.mets,constraintID));
    
    % link metabolic flux and enzyme concentration
    enzID     = 'TPS';  % enzye ID in reaction vector
    model_pa  = addReaction(model_pa,enzID,...
            'reactionName','translationalSectorConcentration',...
            'metaboliteList',{constraintID},...
            'stoichCoeffList',transMolMass*1e-6,...  % unit of enzyme concentration is nmol/g, transform to mass flux mg/g/h
            'reversible',0,...
            'lowerBound',0,...
            'upperBound',1e5,...
            'objectiveCoef',0,...
            'subSystem','Translational Sector',...
            'printLevel',0);
    TPS_pos     = find(strcmp(model_pa.rxns,enzID)); 
    
    % add linear connection to growth rate
    model_pa.S(transSectConstraint_pos,transSector.bmRxnNum) = -transSector.TPS_mu*1000;
       
    % add to total protein constraint
    model_pa.S(TPC_pos,TPS_pos)    = transMolMass*1e-06;
    
    % save info
    opt_pa.translationalSector.molMass      = transMolMass;
    opt_pa.translationalSector.TPS_0        = transSector.TPS_0;
    opt_pa.translationalSector.TPS_mu       = transSector.TPS_mu;
    opt_pa.translationalSector.rxnID        = enzID;
    opt_pa.translationalSector.cnstrID      = constraintID;
    opt_pa.translationalSector.bmRxnID      = transSector.bmRxnID;
end


%% Add excess enzymes protein sector
if ~isempty(excEnzSector)
    fprintf('Add excess enzymes sector ...\n')
    % define molar mass of a concentration unit of the excess enzymes sector
    excEnzMolMass_default    = 3.947778784340140e04; % mean enzymes mass E. coli [g/mol]
    if ~isfield(excEnzSector,'molMass')
        excEnzMolMass    = excEnzMolMass_default;
    else
        if ~isempty(excEnzSector.molMass)
            excEnzMolMass    = excEnzSector.molMass;
        else
            excEnzMolMass    = excEnzMolMass_default;
        end
    end
    
    
    constraintID    = 'EEPS';    % constraint ID (goes as metabolite ID)
    model_pa  = addMetabolite(model_pa,constraintID,...
                'metName','excessEnzymesSector',...
                'b',excEnzSector.EEPS_0*1000,'csense','G');  
    excessEnzymeConstraint_pos     = find(strcmp(model_pa.mets,constraintID));
    
    enzID     = 'EEPS';
    model_pa  = addReaction(model_pa,enzID,...
            'reactionName',enzID,...
            'metaboliteList',{constraintID},...
            'stoichCoeffList',excEnzMolMass*1e-6,...  % unit of enzyme concentration is nmol/g
            'reversible',0,...
            'lowerBound',0,...
            'upperBound',1e6,...
            'objectiveCoef',0,...
            'subSystem','Excess Enzymes Sector',...
            'printLevel',0);
    EEPS_pos    = find(strcmp(model_pa.rxns,enzID)); 
    
    % connect enzyme sector to substrate uptake reaction 
    model_pa.S(excessEnzymeConstraint_pos,excEnzSector.subsRxnNum)    ...
                    = (excEnzSector.EEPS_0*1000)/excEnzSector.subsUptakeMax;
    
    % add to total protein constraint
    model_pa.S(TPC_pos,EEPS_pos)    = excEnzMolMass*1e-06;            
                
    % save info
    opt_pa.excessEnzymesSector.subsRxnID        = {excEnzSector.subsRxnID};
    opt_pa.excessEnzymesSector.subsRxnNum       = excEnzSector.subsRxnNum;
    opt_pa.excessEnzymesSector.subsUptakeMax    = excEnzSector.subsUptakeMax;
    opt_pa.excessEnzymesSector.EEPS_0           = excEnzSector.EEPS_0;
    opt_pa.excessEnzymesSector.molMass          = excEnzMolMass;
    opt_pa.excessEnzymesSector.cnstrID          = constraintID;
    opt_pa.excessEnzymesSector.rxnID            = enzID;
end

%% Add custom protein sectors
if ~isempty(customSectors)
    for i=length(customSectors)
        fprintf(['Add custom protein sector ',num2str(i),' ...\n'])
        customSector    = customSectors{i};     % choose sector
        % define molar mass of a concentration unit of the excess enzymes sector
        cstmMolMass_default     = 3.947778784340140e04; % mean enzymes mass E. coli [g/mol]
        if ~isfield(customSector,'molMass')
            cstmMolMass    = cstmMolMass_default;
        else
            if ~isempty(customSector.molMass)
                cstmMolMass    = customSector.molMass;
            else
                cstmMolMass    = cstmMolMass_default;
            end
        end
        
        % define name of the sector
        if isfield(customSector,'name')
            constraintID    = ['CPS_',customSector.name];
            constraintName  = ['customProteinSector_',customSector.name];
            enzID           = ['CPS_',customSector.name];
            enzName         = ['customProteinSector_',customSector.name];
        else
            constraintID    = ['CPS_',num2str(i)];
            constraintName  = ['customProteinSector_',num2str(i)];
            enzID           = ['CPS_',num2str(i)];
            enzName         = ['customProteinSector_',num2str(i)];
        end
        
        model_pa  = addMetabolite(model_pa,constraintID,...
                    'metName',constraintName,...
                    'b',customSector.CPS_0*1000,'csense','E');  
        customConstraint_pos    = find(strcmp(model_pa.mets,constraintID));

        model_pa  = addReaction(model_pa,enzID,...
                'reactionName',enzName,...
                'metaboliteList',{constraintID},...
                'stoichCoeffList',cstmMolMass*1e-6,...  % unit of enzyme concentration is nmol/g
                'reversible',0,...
                'lowerBound',0,...
                'upperBound',1e6,...
                'objectiveCoef',0,...
                'subSystem','Custom Protein Sector',...
                'printLevel',0);
        CPS_pos     = find(strcmp(model_pa.rxns,enzID)); 
    
        % connect enzyme sector to substrate uptake reaction 
        model_pa.S(customConstraint_pos,customSector.linRxnNum)    ...
                    = -(customSector.CPS_s*1000);
    
        % add to total protein constraint
        model_pa.S(TPC_pos,CPS_pos)    = cstmMolMass*1e-06;    
        
        % save info
        
        opt_pa.(constraintName).linRxnID    = customSector.linRxnID;
        opt_pa.(constraintName).linRxnNum   = customSector.linRxnNum;
        opt_pa.(constraintName).CPS_0       = customSector.CPS_0;
        opt_pa.(constraintName).CPS_s       = customSector.CPS_s;
        opt_pa.(constraintName).molMass     = cstmMolMass;
        opt_pa.(constraintName).cnstrID     = constraintID;
        opt_pa.(constraintName).rxnID       = enzID;
    end 
end



%% finish model
model_pa.opt_pa     = opt_pa;

%% test model
fprintf('Test protein allocation model ...\n')
% setup model
sol         = optimizeCbModel(model_pa,'max');
if strcmp(sol.origStat,'INFEASIBLE')
    warning('Protein allocation model test case is infeasible')
else
    fprintf('Protein allocation model test case was successful.\n')
    fprintf(['Optimal objective function value: ',num2str(sol.f),'\n'])
end


%% Additional, embedded functions

fprintf('Done\n')
end