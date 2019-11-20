%% Example script for builiding and handling proein allocation models
% 
% 
% 
% 
% Author: Tobias B. Alter
% NOV 15, 2019

clear
%% load model
load('Example data\iML1515_irreversible.mat');


%% load data
% active enzyme sector
[~,~,RAW]   = xlsread('Example data\proteinAllocationModel_iML1515_EnzymaticData.xls','ActiveEnzymes');
actEnzSector    = [];

[row,col]   = find(strcmp(RAW,'rxnID'));
actEnzSector.rxnID  = RAW((row+1):end,col);

[row,col]   = find(strcmp(RAW,'kcat'));
actEnzSector.kcat  = [RAW{(row+1):end,col}]';

[row,col]   = find(strcmp(RAW,'molMass'));
actEnzSector.molMass  = [RAW{(row+1):end,col}]';

% excess enzymes sector
[~,~,RAW]   = xlsread('Example data\proteinAllocationModel_iML1515_EnzymaticData.xls','ExcessEnzymes');
excEnzSector    = [];

[row,col]   = find(strcmp(RAW,'subsRxnID'));
excEnzSector.subsRxnID  = RAW{row,col+1};

[row,col]   = find(strcmp(RAW,'subsUptakeMax'));
excEnzSector.subsUptakeMax  = RAW{row,col+1};

[row,col]   = find(strcmp(RAW,'EEPS_0'));
excEnzSector.EEPS_0  = RAW{row,col+1};

[row,col]   = find(strcmp(RAW,'molMass'));
excEnzSector.molMass  = RAW{row,col+1};

% translational protein sector
[~,~,RAW]   = xlsread('Example data\proteinAllocationModel_iML1515_EnzymaticData.xls','Translational');
transSector    = [];

[row,col]   = find(strcmp(RAW,'bmRxnID'));
transSector.bmRxnID  = RAW{row,col+1};

[row,col]   = find(strcmp(RAW,'TPS_0'));
transSector.TPS_0  = RAW{row,col+1};

[row,col]   = find(strcmp(RAW,'TPS_mu'));
transSector.TPS_mu  = RAW{row,col+1};

[row,col]   = find(strcmp(RAW,'molMass'));
transSector.molMass  = RAW{row,col+1};

%% additional parameter
% total condition-dependent protein concentration
totProtConc     = 0.258;

%% build protein allocation model
model_pa = buildPAM(model_i,totProtConc,...
                'actEnzSector',actEnzSector,...
                'excEnzSector',excEnzSector,...
                'transSector',transSector);


%% change parameter values
model_pa    = changePAMParameter(model_pa,...
                'kcat',[1,2],'rxnID',{'PDH','PGI_f'},...    % change kcat values
                'totProtConc',0.311,...     % change total protein concentration
                'subsRxnID','EX_ac_e_b','subsUptakeMax',20,...  % change substrate
                'EEPS_0',excEnzSector.EEPS_0*0.9,...    % change protein allocation towards the excess enzymes sector
                'TPS_0',transSector.TPS_0*1.1,...   % change protein allocation towards the translational sector
                'TPS_mu',transSector.TPS_mu*1.5,... % change slope of the translational protein sector
                'printFlag',1);     % change print flag



%% test model
sol     = optimizeCbModel(model_pa,'max');


