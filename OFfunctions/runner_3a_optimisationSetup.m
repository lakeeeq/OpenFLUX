%%%specify measurement data%%%
%%%modify optimisation parameters%%%
%%%run optimisation%%%%

OFspec = OpenFLUX.OFobjSpecification(mfilename,OFspecFileName);

load(OFspec.modelObjSaveName);

OF.isOptimisation = true;
OF.genLabelledSubstrate;
OF.metDataFileName = OFspec.metDataFileName;

% %%%%quick test by loading existing feasible soln%%%%%%
% if OF.isODEsolver
%     aa = load('ODE_xFeas_op');
% else
%     aa = load('SBR_xFeas_op');
% end
% %%%%%%%%%%%%%%%%%%%%


%%%setup data input%%%
%run this for the first time to import and estimate data error
%to avoid accidental ovewrite, rename saved file to *_keep.mat
%reload same dataset to ensure optimisation is consistently parameterised  (initial conc)
if OFspec.reLoadData
    OF.reEstimateError('load',[]);%%reload previous error estimates
    if isempty(OF.dataMet)
        disp('failed to reload data')
        return
    end
else
    OF.importMetData;
    OF.reEstimateError('generate',10000);%%run MC to estimate error
    OF.reEstimateError('save',[]);
end

%specify extra constraints
OFadditionalConstraintsScript

%%%%setup optimisation%%%%
OF.intKntPos = sort(rand(size(OF.intKntPos)));%%%guess rand knot
% OF.intKntPos = aa.xKnot;%%%reuse solution
simParas = OF.prepSimulation;

%%%%%%%%%%%
%specify extra data (radiation)
OFadditionalDataScript

opSave.lb = simParas.lb;
opSave.ub = simParas.ub;
opSave.xFeas = [];
opSave.solnIndex = 0;
opSave.xFitSeries = struct();%datetime saved, xStart, tElapse, xFinish, fval, exitflag
opSave.metDataFileName = OF.metDataFileName;
opSave.bumpUpXFeasIniConc = OFspec.bumpUpXFeasIniConc;
    
%%%%run this to generate many optimisation instances for HPC
if OFspec.isForHPC
    fileListHPC = {};
    disp('creating and saving HPC optimisation instances');
    noItt = OFspec.HPCitt;
else
    noItt = 1;
end
for i = 1:noItt
    OF.intKntPos = sort(rand(size(OF.intKntPos)));%%%guess rand knot
    simParas = OF.prepSimulation;
    x0 = rand(size(simParas.lb));
    fitFxn = OF.generateFitFxn(additionalData);
    opSave.datetimeCreated = datetime('now','Format','yyyyMMdd_HHmmSSS');
    opSave.fitFxn = fitFxn;
    opSave.x0 = x0;
    opSave.xIntKnot = OF.intKntPos;
    
    if OF.isODEsolver
        opSave.isODE = true;
        opSave.conFxn = simParas.conFxn;
        opSave.saveFileName = strcat(['ODEop_' char(opSave.datetimeCreated) '.mat']);
    else
        opSave.isODE = false;
        opSave.stepBTWsample = OF.stepBTWsample;
        opSave.AconParas = simParas.AconParas;
        opSave.saveFileName = strcat(['SBRop_' char(opSave.datetimeCreated) '.mat']);
    end
    if OFspec.isForHPC
        fileListHPC{i,1} = opSave.saveFileName;
    end
    save(strcat(OFspec.opSaveFolder, opSave.saveFileName), 'opSave','OF');
    disp(strcat(OFspec.opSaveFolder, opSave.saveFileName));
    pause(1);
end
if OFspec.isForHPC
    save HPCopFileList fileListHPC
    disp('please rename HPCopFileList.mat');
else
    disp('optimisation file created');
end
