%%%%change step size of SBR%%%
%%%load and change specific (listed) optimisation instances%%%
clearvars
clc
newSampleSteps = [500	400	200	75	75	50];
opSaveFolder = 'OPinstances/';
fileListToChange  = {
'SBRop_20181116_1038005.mat'
'SBRop_20181116_1038243.mat'
'SBRop_20181116_1038346.mat'
'SBRop_20181116_1038497.mat'
'SBRop_20181116_1038632.mat'
'SBRop_20181116_1038870.mat'
'SBRop_20181116_1039045.mat'
'SBRop_20181116_1039433.mat'
'SBRop_20181116_1039703.mat'
'SBRop_20181116_1039791.mat'
};

load OFobj
% OF.isODEsolver = true;%toggle off to choose SBR
OF.isOptimisation = true;
OF.genLabelledSubstrate;

%constrain accoa_out and glycogen_out initial to near zero 0.01
hitMet = strcmp('ACCOA_out',OF.metListInt);
OF.concScale(hitMet,:) = 0.01;
hitMet = strcmp('GLYCOGEN_out',OF.metListInt);
OF.concScale(hitMet,:) = 0.01;

for i = 1:numel(fileListToChange)
    load(strcat([opSaveFolder fileListToChange{i}]));
    OF.metDataFileName = opSave.metDataFileName;%for insulin
    OF.reEstimateError('load',[]);
    OF.stepBTWsample = newSampleSteps;
    OF.intKntPos = opSave.xIntKnot;
    simParas = OF.prepSimulation;
    
    %%%%%%%%%%%
    %specify extra data (radiation)
    additionalDataScript
    %%%%%%%%%%%
    
    fitFxn = OF.generateFitFxn(additionalData);
    opSave.fitFxn = fitFxn;
    opSave.AconParas = simParas.AconParas;
    opSave.stepBTWsample = OF.stepBTWsample;
    save(strcat(opSaveFolder, opSave.saveFileName), 'opSave');
end
