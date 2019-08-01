OFspec.modelObjSaveName = 'OFobj_SBR';
% OFspec.isODEsolver = true;%choose ODE
OFspec.isODEsolver = false;%choose SBR
switch scriptCalling
    case 'runner_1_modelSetup'
        OFspec.modelCondition = 'insulin';%tag for OF object
        OFspec.modelFileName = 'model.txt';%metabolic model
        OFspec.ionFormFileName = 'metIonFormula.txt';%formula for MZs, natural isotope correction
        OFspec.natEndo13Cenrich = 0.0107;%endogenous metabolites' natural carbon enrichment
        OFspec.natSub13Cenrich = 0.0107;%input substrates' natural carbon enrichment
        OFspec.intKntPos = [0.2];%specify bspline internal knots
        OFspec.orderS = 3;%specify bspline order
        OFspec.labelledSub = {'GLC_in' [0.99 0.99 0.99 0.99 0.99 0.99]};%%%specify positional enrichment
        OFspec.sampleTime = [1 5 10 20 40 60];%indicate sampling time points (still required for simulations)
                
        OFspec.stepBTWsample = [20 20 20 20 20 20];%required for SBR, low res
        % OFspec.stepBTWsample = [500 400 200 75 75 50];%required for SBR, high res
        
        OFspec.odeSimTime = [];%required for ODE15s, sampleTime must be within this range. set as empty if using SBR
        
        %%%optimisation parameters, current showing default
        OFspec.concBound = [0.01 10000];
        OFspec.fluxBound = [10 100000];        
        
        
    case 'runner_2_genFeasible'
        OFspec.simSaveFolder = 'SIMinstances/';%folder name to load/save all the optimisation instances
        
        %         OFspec.bumpUpXFeasIniConc = false;%use this to overcome negatives in ODE solver by increasing error tolerance
        OFspec.bumpUpXFeasIniConc = true;%use this to overcome negatives in ODE solver by increasing initial metabolite concentration
        
        
    case 'runner_3a_optimisationSetup'
        OFspec.opSaveFolder = 'OPinstances/';%folder name to load/save all the optimisation instances
        
        OFspec.metDataFileName = 'metDataIns';%for insulin
        %         OFspec.metDataFileName = 'metDataBas.txt';%for basal
        
        OFspec.reLoadData = true;%re-use generated measurements average and error
        %         OFspec.reLoadData = false;%regenerate measurements average and error
        
        %         OFspec.bumpUpXFeasIniConc = false;%use this to overcome negatives in ODE solver by increasing error tolerance
        OFspec.bumpUpXFeasIniConc = true;%use this to overcome negatives in ODE solver by increasing initial metabolite concentration
        
        %         OFspec.isForHPC = false;%choose this if for local run
        OFspec.isForHPC = true;%choose this if optimisation is performed on computing clusters
        OFspec.HPCitt = 10;%set number of parallel runs on computing clusters
        
        
    case 'runner_3c_changeStepSize'
        OFspec.newSampleSteps = [500 400 200 75 75 50];
        OFspec.opSaveFolder = 'OPinstances/';
        OFspec.fileListToChange  = {
            'ODEop_20190801_1436882'
            'SBRop_20190801_1441670'
            'SBRop_20190801_1441354'
            'SBRop_20190801_1441944'
            };
        
        
    case 'runner_3d_changeSBRtoODE'
        OFspec.isODEsolver = true;%now change SBR to ODE
        %         OFspec.isODEsolver = false;%now change ODE to SBR
        OFspec.opSaveFolder = 'OPinstances/';%folder name to load/save all the optimisation instances
        OFspec.fileListToChange  = {
            'ODEop_20190801_1436882'
            'SBRop_20190801_1441670'
            'SBRop_20190801_1441354'
            'SBRop_20190801_1441944'
            };
        OFspec.odeSimTime = [0:1:60];%required for ODE15s, sampleTime must be within this range.
        
        
    case 'runner_4_visualiseSoln'
%         OFspec.loadOP = false;
%         OFspec.loadFolder = 'SIMinstances/';
%         OFspec.fileName = 'SBRsim_20190801_1559309';
        
        OFspec.loadOP = true;
        OFspec.loadFolder = 'OPinstances/';
        OFspec.fileName = 'SBRop_20190801_1436634';
        OFspec.simResidualError = false;
        
end