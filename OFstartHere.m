clearvars
clc
addpath OFfunctions %can be added to MATLAB startup
addpath inputs %can be added to MATLAB startup
% rng('shuffle');%run this once upon starting MATLAB

% OFspecFileName = 'OFspec_ODE';
OFspecFileName = 'OFspec_SBR';

% runner_1_modelSetup

% runner_2_genFeasible

% runner_3a_optimisationSetup
% runner_3b_optimisationRun('OPinstances/','SBRop_20190801_1436634');

% for i = 1:10
%     runner_3b_optimisationRun_HPC(i,'HPCopFileList_SBR','OPinstances/');
% end

% runner_3c_changeStepSize

% runner_3d_changeSBRtoODE %also can change from ODE to SBR

runner_4_visualiseSoln
return



%%%%%other useful functions%%%%%
%%%unpac variablesk%%%
varName = fieldnames(opInput);
for i = 1:numel(varName)
    eval([varName{i} '=opInput.' varName{i} ';']);
end



%%%get values of private properties%%%
OF.getPrivProp('property name');
%replace "property name" with variable name. For example OF.getPrivProp('bigEMUmodel');
%or leave empty to get list of private properties
