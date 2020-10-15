%%%specify measurement data%%%
%%%modify optimisation parameters%%%
%%%run optimisation%%%%

OFspec = OpenFLUX.OFobjSpecification(mfilename,OFspecFileName);

load(OFspec.modelObjSaveName);

if ~isfolder(OFspec.opSaveFolder)
    mkdir(OFspec.opSaveFolder);
end

OF.isOptimisation = true;
OF.genLabelledSubstrate;
OF.metDataFileName = OFspec.metDataFileName;
OF.inputDirectory = OFspec.inputDirectory;

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
if OF.isDynamic
    OF.intKntPos = sort(rand(size(OF.intKntPos)));%%%guess rand knot
    % OF.intKntPos = aa.xKnot;%%%reuse solution
end

simParas = OF.prepSimulation;

%%%%%%%%%%%
%specify extra data (radiation)
OFadditionalDataScript
OF.additionalData = additionalData;

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
    if OF.isDynamic
        OF.intKntPos = sort(rand(size(OF.intKntPos)));%%%guess rand knot
    end
    simParas = OF.prepSimulation;
    x0 = rand(size(simParas.lb));
    fitFxn = OF.generateFitFxn;
    opSave.datetimeCreated = datetime('now','Format','yyyyMMdd_HHmmSSS');
    opSave.fitFxn = fitFxn;
    opSave.x0 = x0;
    opSave.xIntKnot = OF.intKntPos;
    
    if OF.isODEsolver
        opSave.isODE = true;
        opSave.isDynamic = true;
        opSave.conFxn = simParas.conFxn;
        opSave.saveFileName = strcat(['ODEop_' char(opSave.datetimeCreated) '.mat']);
    elseif OF.isDynamic
        opSave.isODE = false;
        opSave.isDynamic = true;
        opSave.stepBTWsample = OF.stepBTWsample;
        opSave.AconParas = simParas.AconParas;
        opSave.saveFileName = strcat(['SBRop_' char(opSave.datetimeCreated) '.mat']);
    else
        opSave.isODE = false;
        opSave.isDynamic = false;
        opSave.conFxn = simParas.conFxn;
        opSave.saveFileName = strcat(['SSop_' char(opSave.datetimeCreated) '.mat']);
    end
    if OFspec.isForHPC
        fileListHPC{i,1} = opSave.saveFileName;
    end
    save(strcat(OFspec.opSaveFolder, opSave.saveFileName), 'opSave','OF');
    disp(strcat(OFspec.opSaveFolder, opSave.saveFileName));
    pause(0.2);
end
if OFspec.isForHPC
    save HPCopFileList fileListHPC
    disp('please rename HPCopFileList.mat');
else
    disp('optimisation file created');
end

disp('list of internal metabolites:')
for i = 1:size(OF.metListInt)
    fprintf('%1.0f\t%s\n',i,OF.metListInt{i});    
end
fprintf('\n');

fprintf('list of metabolite data:\n')
fprintf('EMUname(col2)\tmetRow(col3)\t[EMUindexMap](col4)\n');
for i = 1:size(OF.dataMet)
    fprintf('%1.0f\t%s\t%1.0f\t[%1.0f,%1.0f]\n',i,OF.dataMet{i,1},OF.dataMet{i,8},OF.dataMet{i,10}(1),OF.dataMet{i,10}(2));    
end
fprintf('\n');
