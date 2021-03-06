addpath(pwd);
pwdString = pwd;
hitFS = regexp(pwdString,filesep);
OFspec.inputDirectory = pwdString((hitFS(end)+1):end);

%%
switch scriptCalling
    case 'runner_1_modelSetup'
        OFspec.modelObjSaveName = 'OFobj_ODEdynaToy';
        OFspec.isODEsolver = true;%choose ODE
        %         OFspec.isODEsolver = false;%choose SBR
        OFspec.isDynamic = true;
        
        OFspec.modelCondition = 'dynamicToyODE';%tag for OF object
        OFspec.modelFileName = 'model.txt';%metabolic model
        %         OFspec.ionFormFileName = 'metIonFormula.txt';%formula for MZs, natural isotope correction
        OFspec.intKntPos = [0.1 0.3 0.5 0.6];%specify bspline internal knots
        OFspec.orderS = 3;%specify bspline order
        OFspec.labelledSub = {'PYR_in' [0.99 0.99 0.99] 1};%%%specify positional enrichment
        OFspec.sampleTime = [1 5 10 15 20 25 30 40 50 60];%indicate sampling time points (still required for simulations)
        
        %         OFspec.stepBTWsample = [20 20 20 20 20 20 20 20 20 20];%required for SBR, low res
        
        OFspec.odeSimTime = [0:1:60];%required for ODE15s, sampleTime must be within this range. set as empty if using SBR
        
        %%%additional parameters, currently showing default
        OFspec.concBound = [50 200];
        OFspec.fluxBound = [0.01 25];
        OFspec.natEndo13Cenrich = 0.0107;%endogenous metabolites' natural carbon enrichment
        OFspec.natSub13Cenrich = 0.0107;%input substrates' natural carbon enrichment
        
        
        %%
    case 'runner_2_genFeasible'
        OFspec.modelObjSaveName = 'OFobj_ODEdynaToy';
        OFspec.simSaveFolder = 'SIMinstances/';%folder name to load/save all the optimisation instances
        
        %         OFspec.bumpUpXFeasIniConc = false;%use this to overcome negatives in ODE solver by increasing error tolerance
        OFspec.bumpUpXFeasIniConc = true;%use this to overcome negatives in ODE solver by increasing initial metabolite concentration
        OFspec.simulateSoln = true;%run simulation of generated solution, writen into simOutput
        
        
        %%
    case 'runner_3a_optimisationSetup'
        OFspec.modelObjSaveName = 'OFobj_ODEdynaToy';
        OFspec.opSaveFolder = 'OPinstances/';%folder name to load/save all the optimisation instances
        
        OFspec.metDataFileName = 'modelData_toy';%for insulin
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
        OFspec.newSampleSteps = [20 20 20 20 20 20 20 20 20 20];%required for SBR, low res
        
        OFspec.opSaveFolder = 'OPinstances/';%folder name to load/save all the optimisation instances
        OFspec.fileListToChange  = {
            'SBRop_20190801_1441942'
            'ODEop_20190801_1436882'
            };
        
        
        %%
    case 'runner_4a_visualiseSoln'
        %         OFspec.loadFolder = 'SIMinstances/';
        %         OFspec.fileName = 'ODEsim_20200805_1328042';
        
        OFspec.loadFolder = 'OPinstances/';
        OFspec.fileName = 'ODEop_20200805_1332751';
        OFspec.simResidualError = true;
        
        
        %%
    case 'runner_4b_visualiseMetData'
        OFspec.loadFolder = 'OPinstances/';
        OFspec.fileName = 'ODEop_20200805_1332751';
        
        
        %%
    case 'runner_5a_MonteCarloSetup'
        OFspec.loadFolder = 'OPinstances/';
        OFspec.fileName = 'ODEop_20190802_1349899';%OP instance to clone
        OFspec.noCase = 3;%number of data corruption
        OFspec.noRepeatsPerCase = 5;%how many repeated (multistart) optimisation for each corruption
        OFspec.mcSavFolder = 'MCinstances/';
        
        
        %%
    case 'runner_6_compileVisualiseOps'
        OFspec.loadFolder = 'MCinstances/';
        OFspec.opReferenceFile = 'SBRmc_20190903_1324165';
        OFspec.compileFileName = 'compiledOpResults_1';
        OFspec.reloadCompiledResults = false;
        OFspec.mcBestSoln = false;
        %%
end