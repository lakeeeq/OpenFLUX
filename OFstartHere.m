clearvars
clc
addpath OFfunctions %can be added to MATLAB startup
% rng('shuffle');%run this once upon starting MATLAB

%%%to do next, features to add%%%
%compile optimised solns from OP/MC runs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% object specification
% OFspecFileName = 'inputs_DYNAMICtoy/OFspec_ODEdynaToy';
% OFspecFileName = 'inputs_DYNAMICtoy/OFspec_SBRdynaToy';
% OFspecFileName = 'inputs_adipocytes/OFspec_ODEadipocyte';
% OFspecFileName = 'inputs_adipocytes/OFspec_SBRadipocyte';
OFspecFileName = 'inputs_SStoy/OFspec_SStoy';

%% task to run
taskToDO = 'runner_1';%model build and setup (for SS too)
% taskToDO = 'runner_2';%generate feasible solutions and simulate (for SS too)
% taskToDO = 'runner_3a';%optimisation setup (for SS too)
% taskToDO = 'runner_3c';%change SBR step size
% taskToDO = 'runner_3d';%change SBR to ODE, and vice-versa
% taskToDO = 'runner_4a';%visualise flux solution
% taskToDO = 'runner_4b';%visualise metabolite data
% taskToDO = 'runner_5';%set up monte carlo runs(MC instances, for SS too)

%{
taskToDO = 'runner_3b';%run optimisation (single instance, for SS too)
OFspec.saveFolder = 'OPinstances/';
OFspec.OPinstance = 'ODEop_20190903_1321347';
%}

%{
taskToDO = 'runner_3b_hpc';%run optimisations on HPC (for testing, for MC, for SS too)
OPspec.noItt = 50;
OFspec.saveFolder = 'MCinstances/';
OFspec.fileListVarSaved = 'HPCmcFileList_1';
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


switch taskToDO
    case 'runner_1'
        runner_1_modelSetup
    case 'runner_2'
        runner_2_genFeasible
    case 'runner_3a'
        runner_3a_optimisationSetup
    case 'runner_3b'
        runner_3b_optimisationRun(OFspec.saveFolder,OFspec.OPinstance);
    case 'runner_3b_hpc'
        for i = 1:OPspec.noItt
            runner_3b_optimisationRun_HPC(i,OFspec.fileListVarSaved,OFspec.saveFolder);
        end
    case 'runner_3c'
        runner_3c_changeStepSize
    case 'runner_3d'
        runner_3d_changeSBRtoODE
    case 'runner_4a'
        runner_4a_visualiseSoln
    case 'runner_4b'
        runner_4b_visualiseMetData
    case 'runner_5'
        runner_5_MonteCarloSetup
end
return


%%%%%other useful functions%%%%%
%read and write into OpenFLUX private properies
OF.getPrivProp('property name');%read
OF.writePrivProp('property name', 'input variable');%write
% replace "property name" with variable name. For example
% newJacOut = OF.getPrivProp('jacOut');
% OF.writePrivProp('jacOut',newJacOut);
EMUmodelOutput = OF.getPrivProp('EMUmodelOutput');
