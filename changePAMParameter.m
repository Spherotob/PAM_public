function model_pa = changePAMParameter(model_pa,varargin)
% Change proein allocation model parameters
% 
% 
% INPUTS:
%     model_pa:     protein allocation model
% 
% 
%     varargin      Optional Inputs provided as 'ParameterName', Value pairs.
% 
%               *kcat       vector containing kcat values to be changed.
%                           Also requires "rxnID" as an input
%                           argument.(unit: 1/s)
% 
%               *rxnID      Character vector of Reaction ID of enzymatic reactions for which
%                           kcat values are to be changed. May be a cell array 
%                           if multiple reactions should be changed. Also requires
%                           "kcat" as an input argument.
% 
%               *totProtConc    New concentration of total
%                               condition-dependent protein (unit: g/g_CDW)
% 
%               *subsRxnID  Character vector with new substrate uptake
%                           reaction ID for parametrization of the excess enzymes sector.
%                           Can also be a cell array if
%                           multiple substrates should be considered at the
%                           same time. Also requires "subsUptakeMax" as an
%                           input argument.
% 
%               *subsUptakeMax      Vecor containing the maximal substrate
%                                   uptake rates for the new substrate
%                                   upake reactions.
% 
%               *EEPS_0     Protein allocated to excess enzymes sector at
%                           zeros substrate uptake (unit: g/g_CDW).
% 
%               *TPS_0      Protein allocated to translational protein
%                           sector at zero growth (unit: g/g_CDW).
% 
%               *TPS_mu     Slope protein allocated to translational sector
%                           per unit increase of growth rate (unit: g h/g_CDW)
% 
% 
% 
% 
%               *printFlag  (1): prints ouput (default); (0) no printed output
% 
% 
% OUTPUT
%   model_pa:   Protein allocation model including parameter changes
%
% 
% Author: Tobias B. Alter
% NOV 15, 2019

%% setup parameter change
% preallocate variables
kcat                = [];
rxnID               = {};
totProtConc         = [];
subsRxnID           = {};
subsUptakeMax       = [];
EEPS_0_parameter    = [];
TPS_0_parameter     = [];
TPS_mu_parameter    = [];

if any(strcmp(varargin,'kcat'))
    kcat    = varargin{find(strcmp(varargin,'kcat'))+1};
end
if any(strcmp(varargin,'rxnID'))
    rxnID    = varargin{find(strcmp(varargin,'rxnID'))+1};
end
if any(strcmp(varargin,'totProtConc'))
    totProtConc    = varargin{find(strcmp(varargin,'totProtConc'))+1};
end
if any(strcmp(varargin,'subsRxnID'))
    subsRxnID    = varargin{find(strcmp(varargin,'subsRxnID'))+1};
end
if any(strcmp(varargin,'subsUptakeMax'))
    subsUptakeMax    = varargin{find(strcmp(varargin,'subsUptakeMax'))+1};
end
if any(strcmp(varargin,'EEPS_0'))
    EEPS_0_parameter    = varargin{find(strcmp(varargin,'EEPS_0'))+1};
end
if any(strcmp(varargin,'TPS_0'))
    TPS_0_parameter    = varargin{find(strcmp(varargin,'TPS_0'))+1};
end
if any(strcmp(varargin,'TPS_mu'))
    TPS_mu_parameter    = varargin{find(strcmp(varargin,'TPS_mu'))+1};
end
if any(strcmp(varargin,'printFlag'))
    printFlag    = varargin{find(strcmp(varargin,'printFlag'))+1};
else
    printFlag   = 1;
end

%% load model data
opt     = model_pa.opt_pa;


%% change kcat values
if ~isempty(kcat) && ~isempty(rxnID)
    if ischar(rxnID)
        rxnID   = {rxnID};
    end
    % check provided data
    if length(kcat)~=length(rxnID)
        warning('Number of provided kcat values does not match the number of reaction IDs')
    else
        % load parameter
        rxns2ECrxns         = opt.activeEnzymesSector.rxns2ECrxns;
        cnstrIDs2ECrxns     = opt.activeEnzymesSector.cnstrIDs2ECrxns;
        ECrxns              = opt.activeEnzymesSector.ECrxns;
        for i=1:length(kcat)
            pos     = find(strcmp(rxns2ECrxns,rxnID{i}));
            model_pa.S(strcmp(model_pa.mets,cnstrIDs2ECrxns{pos}),...
                        strcmp(model_pa.rxns,ECrxns{pos})) ...
                = -(kcat(i)*3600)*1e-6;
            % save new kcat value
            oldkcat                 = opt.activeEnzymesSector.kcat2ECrxns(pos);
            opt.activeEnzymesSector.kcat2ECrxns(pos)    = kcat(i);
            % print
            if printFlag
                fprintf([rxnID{i},': kcat value changed from: ',num2str(oldkcat),' 1/s to ',...
                    num2str(kcat(i)),' 1/s\n'])
            end
        end
    end

end


%% change concentration of total condition-dependent protein
if ~isempty(totProtConc)
    % change parameter
    model_pa.b(strcmp(model_pa.mets,opt.totalProtein.cnstrID))    ...
                = totProtConc*1000;
    % save new value
    oldconc     = opt.totalProteinConcentration;
    opt.totalProtein.totalProteinConcentration   = totProtConc;
    % print
    if printFlag
        fprintf(['Total condition-dependent protein concentration changed from: ',...
            num2str(oldconc),' g/g_CDW to ', num2str(totProtConc),' g/g_CDW \n'])
    end
end

%% change concentration of maximal excess enzyme protein allocation
if ~isempty(EEPS_0_parameter)
    % change parameter
    model_pa.b(strcmp(model_pa.mets,opt.excessEnzymesSector.cnstrID)) ...
                = EEPS_0_parameter*1000;
    % save new value
    oldconc     = opt.excessEnzymesSector.EEPS_0;
    opt.excessEnzymesSector.EEPS_0   = EEPS_0_parameter;
    % print
    if printFlag
        fprintf(['Maximal protein amount allocated to excess enzymes sector changed from: ',...
            num2str(oldconc),' g/g_CDW to ', num2str(EEPS_0_parameter),' g/g_CDW \n'])
    end
end

%% change parameters of the excess enzymes sector
if ~isempty(subsRxnID) && ~isempty(subsUptakeMax)
    % check if single or multiple substrates were supplied
    if ischar(subsRxnID)
        subsRxnID   = {subsRxnID};
    end
    % check if excess enzymes sector exists
    if ~isfield(opt,'excessEnzymesSector')
        warning('There is no excess enzymes sector in the model')
    elseif length(subsRxnID)~=length(subsUptakeMax)
        % check if provided data is consistent
        warning('Number of substrate uptake reaction IDs does not mach the number of substrate uptake rate values.')
    else
        % load parameter
        EEPS_pos    = find(strcmp(model_pa.mets,opt.excessEnzymesSector.cnstrID));
        EEPS_0      = opt.excessEnzymesSector.EEPS_0;
            
        % delete link to current substrate
        for i=1:length(opt.excessEnzymesSector.subsRxnID)
            model_pa.S(EEPS_pos,find(strcmp(model_pa.rxns,...
                opt.excessEnzymesSector.subsRxnID{i})))   = 0;            
        end
            
        
        % add new substrates
        subsRxnNum      = [];
        subsRxnID_s     = {};
        subsUptakeMax_s = [];
        for i=1:length(subsRxnID)
            if ~any(strcmp(model_pa.rxns,subsRxnID{i}))
                warning(['Substrate uptake reaction ',subsRxnID{i},' not found in the model.'])
            else
                subsRxnNum(end+1,1)         = find(strcmp(model_pa.rxns,subsRxnID{i}));
                subsRxnID_s{end+1,1}        = subsRxnID{i};
                subsUptakeMax_s(end+1,1)    = subsUptakeMax(i);
                % change parameter
                model_pa.S(EEPS_pos,subsRxnNum(end)) = (1000*EEPS_0)/subsUptakeMax(i);   
                % print
                if printFlag
                    fprintf(['New substrate uptake reaction: ',subsRxnID{i},'\n'])
                end
            end       
        end
        % save parameter
        opt.excessEnzymesSector.subsRxnNum  = subsRxnNum;
        opt.excessEnzymesSector.subsRxnID   = subsRxnID_s;
        opt.excessEnzymesSector.subsUptakeMax   = subsUptakeMax_s;
    end  
end

%% change protein concentration allocated to the translational sector at zero growth
if ~isempty(TPS_0_parameter)
    % change parameter
    model_pa.b(strcmp(model_pa.mets,opt.translationalSector.cnstrID)) ...    
            = TPS_0_parameter*1000;
    % save new value
    oldconc     = opt.translationalSector.TPS_0;
    opt.translationalSector.TPS_0   = TPS_0_parameter;
    % print
    if printFlag
        fprintf(['Protein allocated to translational sector at zero growth changed from: ',...
            num2str(oldconc),' g/g_CDW to ', num2str(TPS_0_parameter),' g/g_CDW \n'])
    end
end

%% change slope of protein allocated to the translational sector
if ~isempty(TPS_mu_parameter)
    % change parameter
    model_pa.S(strcmp(model_pa.mets,opt.translationalSector.cnstrID), ...  
                strcmp(model_pa.rxns,opt.translationalSector.bmRxnID)) ...
            = -TPS_mu_parameter*1000;
    % save new value
    oldconc     = opt.translationalSector.TPS_mu;
    opt.translationalSector.TPS_mu   = TPS_mu_parameter;
    % print
    if printFlag
        fprintf(['Protein allocation slope of translational sector changed from: ',...
            num2str(oldconc),' g h/g_CDW to ', num2str(TPS_mu_parameter),' g h/g_CDW \n'])
    end
end

% save new options file
model_pa.opt_pa     = opt;
end