%%%specify measurement data%%%
%%%modify optimisation parameters%%%
%%%run optimisation%%%%
clearvars
clc
% rng('shuffle');
opSaveFolder = 'OPinstances/';
load OFobj
OF.isODEsolver = true;%toggle off to choose SBR
OF.isOptimisation = true;
OF.genLabelledSubstrate;
OF.metDataFileName = 'metDataIns';%for insulin
% OF.metDataFileName = 'metDataBas.txt';%for basal

%%%%quick test%%%%%%
if OF.isODEsolver
    aa = load('ODE_xFeas_op');
else
    aa = load('SBR_xFeas_op');
end
%%%%%%%%%%%%%%%%%%%%

reLoadData = true;
% reLoadData = false;
%%%setup data input%%%
%run this for the first time to import and estimate data error
%to avoid accidental ovewrite, rename saved file to *_keep.mat
%reload same dataset to ensure optimisation is consistently parameterised  (initial conc)
if reLoadData    
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

%constrain accoa_out and glycogen_out initial to near zero 0.01
hitMet = strcmp('ACCOA_out',OF.metListInt);
OF.concScale(hitMet,:) = 0.01;
hitMet = strcmp('GLYCOGEN_out',OF.metListInt);
OF.concScale(hitMet,:) = 0.01;

%%%%setup optimisation%%%%

OF.intKntPos = sort(rand(size(OF.intKntPos)));%%%guess rand knot
OF.intKntPos = aa.xKnot;%%%reuse solution
simParas = OF.prepSimulation;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%specify extra data (radiation)
%%%radiolabel data%%%
%%%setting up correction matrix
%this is for FA accumulation data, as acetyl-CoA
%c1 glucose, c2 mix, c3 unlabeled
additionalData{1,1} = pinv([OpenFLUX.cVectGen2([0.01 0.99],2,3)'
    cauchy(OpenFLUX.cVectGen2([1-0.0107 0.0107],1,2)',OpenFLUX.cVectGen2([0.01 0.99],1,2)')
    OpenFLUX.cVectGen2([1-0.0107 0.0107],2,3)']');
%this is for glycogen, as glucose monomer
additionalData{2,1} = pinv([OpenFLUX.cVectGen2([0.01 0.99],6,7) OpenFLUX.cVectGen2([0.9893 0.0107],6,7)]);

%%%specifying data and data position
%%FA accumulation pCmol/mgP accumulated for 1 hr
if strcmp('metDataIns',OF.metDataFileName)%for insulin
    additionalData{1,2} = [76973.028 9680.29035];%average ,se for FA
    additionalData{2,2} = [241.2611 10.16552]*6000;%average ,se for glycogen
else %for basal
    additionalData{1,2} = [16819.87	652.2888155];%average ,se for FA
    additionalData{2,2} = [6.413856	0.478499]*6000;%average ,se for glycogen
end
additionalData{1,3} = OF.findEMUindex('ACCOA_out',[1 1]);
additionalData{2,3} = OF.findEMUindex('GLYCOGEN_out',[1 1 1 1 1 1]);

%%%%run this to generate many optimisation instances for HPC
% %{
opSave.lb = simParas.lb;
opSave.ub = simParas.ub;
opSave.xFeas = [];
opSave.solnIndex = 0;
opSave.xFitSeries = struct();%datetime saved, xStart, tElapse, xFinish, fval, exitflag
for i = 1:10
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
        opSave.AconParas = simParas.AconParas;
        opSave.saveFileName = strcat(['SBRop_' char(opSave.datetimeCreated) '.mat']);
    end
    save(strcat(opSaveFolder, opSave.saveFileName), 'opSave');
end
fileList_HPC = dir('OPinstances\*.mat');%create a file list for deployment, can trim
save fileList_HPC fileList_HPC
return
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = rand(size(simParas.lb));
x0 = aa.xFeas;%%reuse solution
xFeas = x0;

%%%%save optimisation instance%%%%
fitFxn = OF.generateFitFxn(additionalData);

opSave.fitFxn = fitFxn;
opSave.x0 = x0;
opSave.xIntKnot = OF.intKntPos;
opSave.lb = simParas.lb;
opSave.ub = simParas.ub;
opSave.xFeas = [];
opSave.solnIndex = 0;
opSave.xFitSeries = struct();%datetime saved, xStart, tElapse, xFinish, fval, exitflag
opSave.datetimeCreated = datetime('now','Format','yyyyMMdd_HHmmSSS');
if OF.isODEsolver
    opSave.isODE = true;
    opSave.conFxn = simParas.conFxn;
    opSave.saveFileName = strcat(['ODEop_' char(opSave.datetimeCreated) '.mat']);
else
    opSave.isODE = false;
    opSave.AconParas = simParas.AconParas;
    opSave.saveFileName = strcat(['SBRop_' char(opSave.datetimeCreated) '.mat']);
end
save(strcat(opSaveFolder, opSave.saveFileName), 'opSave');
% return%%%terminate here for HPC/batch processing%%%
%%%%continue to test optimisation%%%%
opOptions = optimoptions('fmincon', 'Display','iter','MaxFunEvals',10000);

xGuess = x0;
if opSave.isODE
    while 1
        xFeas = fmincon(@(x)minDistX0(x,opSave.x0),xGuess,[],[],[],[],opSave.lb,opSave.ub,opSave.conFxn,opOptions);
        if all(opSave.conFxn(xFeas)<=0)
            break
        else
            xGuess = xFeas;
        end
    end
else
    [Acon, Bcon] = OpenFLUX.buildSBRcon(opSave.AconParas);
    while 1
        xFeas = fmincon(@(x)minDistX0(x,opSave.x0),xGuess,Acon,Bcon,[],[],opSave.lb,opSave.ub,[],opOptions);
        if all(Acon*xFeas<=Bcon)
            break
        else
            xGuess = xFeas;
        end
    end
end
opSave.xFeas = xFeas;
save(strcat(opSaveFolder, opSave.saveFileName), 'opSave');


%%%running optimisation%%%
disp(fitFxn(xFeas));
return
xStart = xFeas;
tStart = tic;
if opSave.isODE
    [xFinish,fval,exitflag]= fmincon(fitFxn,xStart,[],[],[],[],opSave.lb,opSave.ub,opSave.conFxn,opOptions);
else
    [xFinish,fval,exitflag]= fmincon(fitFxn,xStart,Acon,Bcon,[],[],opSave.lb,opSave.ub,[],opOptions);
end
% xFinish = xStart;fval = 1e11;exitflag = -1;
tElapse = toc(tStart);
opSave.solnIndex = opSave.solnIndex+1;
opSave.xFitSeries(opSave.solnIndex).saveTime = datetime;
opSave.xFitSeries(opSave.solnIndex).xStart = xStart;
opSave.xFitSeries(opSave.solnIndex).xFinish = xFinish;
opSave.xFitSeries(opSave.solnIndex).tElapse = tElapse;
opSave.xFitSeries(opSave.solnIndex).fval = fval;
opSave.xFitSeries(opSave.solnIndex).exitflag = exitflag;
save(strcat(opSaveFolder, opSave.saveFileName), 'opSave');