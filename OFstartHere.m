clearvars
clc
addpath OFfunctions %can be added to MATLAB startup
% rng('shuffle');%run this once upon starting MATLAB

% OFspecFileName = 'OFspec_ODE';
OFspecFileName = 'OFspec_SBR';

% runner_1_modelSetup


% runner_2_genFeasible


% runner_3a_optimisationSetup
% runner_3b_optimisationRun('OPinstances/','SBRop_20190802_1332337');%run inidividual instances
% for i = 1:3
%     runner_3b_optimisationRun_HPC(i,'HPCopFileList_SBR','OPinstances/');
% end


% runner_3c_changeStepSize %for SBR


% runner_3d_changeSBRtoODE %also can change from ODE to SBR


% runner_4a_visualiseSoln
runner_4b_visualiseMetData

% runner_5_MonteCarloSetup
% for i = 1:1
%     runner_3b_optimisationRun_HPC(i,'HPCmcFileList','MCinstances/');
% end

%%%to do next%%%
%compile solns from OP runs
%compile solns from MC runs
return

%%%%%other useful functions%%%%%
%read and write into OpenFLUX private properies
OF.getPrivProp('property name');%read
OF.writePrivProp('property name', 'input variable');%write
%replace "property name" with variable name. For example
%newJacOut = OF.getPrivProp('jacOut');
%OF.writePrivProp('jacOut',newJacOut);




