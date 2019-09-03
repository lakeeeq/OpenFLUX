addpath(pwd);

%%
switch scriptCalling
    case 'runner_1_modelSetup'
        OFspec.modelObjSaveName = 'OFobj_ODE';
        OFspec.isODEsolver = true;%choose ODE
        %         OFspec.isODEsolver = false;%choose SBR
        OFspec.isDynamic = true;
        
        OFspec.modelCondition = 'insulin';%tag for OF object
        OFspec.modelFileName = 'model.txt';%metabolic model
%         OFspec.ionFormFileName = 'metIonFormula.txt';%formula for MZs, natural isotope correction
        OFspec.intKntPos = [0.2];%specify bspline internal knots
        OFspec.orderS = 3;%specify bspline order
        OFspec.labelledSub = {'GLC_in' [0.99 0.99 0.99 0.99 0.99 0.99] 1};%%%specify positional enrichment
        OFspec.sampleTime = [1 5 10 20 40 60];%indicate sampling time points (still required for simulations)
        
        %         OFspec.stepBTWsample = [20 20 20 20 20 20];%required for SBR, low res
        % OFspec.stepBTWsample = [500 400 200 75 75 50];%required for SBR, high res
        
        OFspec.odeSimTime = [0:1:60];%required for ODE15s, sampleTime must be within this range. set as empty if using SBR
        
        %%%additional parameters, currently showing default
        OFspec.concBound = [0.01 10000];
        OFspec.fluxBound = [10 100000];
        OFspec.natEndo13Cenrich = 0.0107;%endogenous metabolites' natural carbon enrichment
        OFspec.natSub13Cenrich = 0.0107;%input substrates' natural carbon enrichment
   
        
 %%       
    case 'runner_2_genFeasible'
        OFspec.modelObjSaveName = 'OFobj_ODE';
        OFspec.simSaveFolder = 'SIMinstances/';%folder name to load/save all the optimisation instances
        
        %         OFspec.bumpUpXFeasIniConc = false;%use this to overcome negatives in ODE solver by increasing error tolerance
        OFspec.bumpUpXFeasIniConc = true;%use this to overcome negatives in ODE solver by increasing initial metabolite concentration
        OFspec.simulateSoln = true;%run simulation of generated solution, writen into simOutput
      
        
%%        
    case 'runner_3a_optimisationSetup'
        OFspec.modelObjSaveName = 'OFobj_ODE';
        OFspec.opSaveFolder = 'OPinstances/';%folder name to load/save all the optimisation instances
        
        OFspec.metDataFileName = 'metDataIns';%for insulin
        %         OFspec.metDataFileName = 'metDataBas.txt';%for basal
        
%         OFspec.reLoadData = true;%re-use generated measurements average and error
                OFspec.reLoadData = false;%regenerate measurements average and error
        
        %         OFspec.bumpUpXFeasIniConc = false;%use this to overcome negatives in ODE solver by increasing error tolerance
        OFspec.bumpUpXFeasIniConc = true;%use this to overcome negatives in ODE solver by increasing initial metabolite concentration
        
        OFspec.isForHPC = false;%choose this if for local run
        %         OFspec.isForHPC = true;%choose this if optimisation is performed on computing clusters
        %         OFspec.HPCitt = 5;%set number of parallel runs on computing clusters
        
        
%%        
    case 'runner_3c_changeStepSize' %%only applies to SBR
        OFspec.newSampleSteps = [500 400 200 75 75 50];
        OFspec.opSaveFolder = 'OPinstances/';
        OFspec.fileListToChange  = {
            'SBRop_20190802_1332437'
            'SBRop_20190802_1249870'
            };
        
        
%%        
    case 'runner_3d_changeSBRtoODE'
        %         OFspec.isODEsolver = true;%now change SBR to ODE
        %         OFspec.odeSimTime = [0:1:60];%required for ODE15s, sampleTime must be within this range.
        
        OFspec.isODEsolver = false;%now change ODE to SBR
        OFspec.newSampleSteps = [20 20 20 20 20 20];
        
        OFspec.opSaveFolder = 'OPinstances/';%folder name to load/save all the optimisation instances
        OFspec.fileListToChange  = {
            'SBRop_20190801_1441942'
            'ODEop_20190801_1436882'
            };
        

%%        
    case 'runner_4a_visualiseSoln'
        OFspec.loadFolder = 'SIMinstances/';
        OFspec.fileName = 'ODEsim_20190801_1654528_keep';
        
%         OFspec.loadFolder = 'OPinstances/';
%         OFspec.fileName = 'SBRop_20190801_1436634';
%         OFspec.simResidualError = false;
        

%%
    case 'runner_4b_visualiseMetData'
        OFspec.loadFolder = 'OPinstances/';
        OFspec.fileName = 'ODEop_20190802_1349899';
            
          
%%      
    case 'runner_5_MonteCarloSetup'
        OFspec.loadFolder = 'OPinstances/';
        OFspec.fileName = 'ODEop_20190802_1349899';%OP instance to clone
        OFspec.noCase = 3;%number of data corruption
        OFspec.noRepeatsPerCase = 5;%how many repeated (multistart) optimisation for each corruption
        OFspec.mcSavFolder = 'MCinstances/';
%%        
end
