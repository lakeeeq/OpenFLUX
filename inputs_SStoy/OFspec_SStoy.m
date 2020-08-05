addpath(pwd);
%%
switch scriptCalling
    case 'runner_1_modelSetup'
        OFspec.modelObjSaveName = 'OFobj_SS';
        OFspec.isDynamic = false;
        
        OFspec.modelCondition = '2009toy';%tag for OF object
        OFspec.modelFileName = 'modelInput_sample.txt';%metabolic model
%         OFspec.ionFormFileName = 'metIonFormula.txt';%formula for MZs, natural isotope correction
        
        OFspec.labelledSub = {%'PYR_EX' [0.99 0.0107 0.0107] 0.5
            'PYR_EX' [0.99 0.99 0.99] 1
            'GLU_EX' [0.99 0.0107 0.0107 0.0107 0.0107] 1};%%%specify positional enrichment
        
        OFspec.par = 100;%weight for hyperfunction
        OFspec.fluxBound = [0.001 1000];
        OFspec.natEndo13Cenrich = 0.0107;%endogenous metabolites' natural carbon enrichment
        OFspec.natSub13Cenrich = 0.0107;%input substrates' natural carbon enrichment
        OFspec.unlabelledFraction = {'LYSX#111111'};%ISA, for endogenous naturally enrichment. otherwise leave empty
        OFspec.unlabelledFraction = {};
   
        
%%       
    case 'runner_2_genFeasible'
        OFspec.modelObjSaveName = 'OFobj_SS';
        OFspec.simSaveFolder = 'SIMinstances/';%folder name to load/save all the optimisation instances
        OFspec.simulateSoln = true;%run simulation of generated solution, writen into simOutput
        
        
%%        
    case 'runner_3a_optimisationSetup'
        OFspec.modelObjSaveName = 'OFobj_SS';
        OFspec.opSaveFolder = 'OPinstances/';%folder name to load/save all the optimisation instances
        
        OFspec.metDataFileName = 'modelData_toy';%data for toy model
        
        %         OFspec.reLoadData = true;%re-use generated measurements average and error
        OFspec.reLoadData = false;%regenerate measurements average and error
        
        OFspec.bumpUpXFeasIniConc = false;%use this to overcome negatives in ODE solver by increasing error tolerance
        %         OFspec.bumpUpXFeasIniConc = true;%use this to overcome negatives in ODE solver by increasing initial metabolite concentration
        
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
        OFspec.isODEsolver = true;%now change SBR to ODE
        OFspec.odeSimTime = [0:1:60];%required for ODE15s, sampleTime must be within this range.
        
        %         OFspec.isODEsolver = false;%now change ODE to SBR
        %         OFspec.newSampleSteps = [20 20 20 20 20 20];
        
        OFspec.opSaveFolder = 'OPinstances/';%folder name to load/save all the optimisation instances
        OFspec.fileListToChange  = {
            'ODEop_20190801_2311205'
            'SBRop_20190802_1332437'
            };
        
 
%%        
    case 'runner_4a_visualiseSoln'
        %         OFspec.loadFolder = 'SIMinstances/';
        %         OFspec.fileName = 'SBRsim_20190801_1559309';
        
        OFspec.loadFolder = 'OPinstances/';
        OFspec.fileName = 'SBRop_20190802_1332337';
        OFspec.simResidualError = true;
        
        
%%       
    case 'runner_4b_visualiseMetData'
        OFspec.loadFolder = 'OPinstances/';
        OFspec.fileName = 'SBRop_20190802_1332337';
        
        
%%        
    case 'runner_5a_MonteCarloSetup'
        OFspec.loadFolder = 'OPinstances/';
        OFspec.fileName = 'SSop_20190830_2302174';%OP instance to clone
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