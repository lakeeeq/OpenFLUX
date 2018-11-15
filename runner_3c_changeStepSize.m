%%%%change step size of SBR%%%
%%%load and change specific (listed) optimisation instances%%%
clearvars
clc
newSampleSteps = [500	400	200	75	75	50];
opSaveFolder = 'OPinstances/';
fileListToChange  = {
%     'SBRop_20181116_0036120.mat'
    'SBRop_20181116_0036182.mat'
    'SBRop_20181116_0036291.mat'
%     'SBRop_20181116_0036448.mat'
    'SBRop_20181116_0036572.mat'
    'SBRop_20181116_0036619.mat'
    'SBRop_20181116_0036823.mat'
    'SBRop_20181116_0036931.mat'
%     'SBRop_20181116_0036993.mat'
    'SBRop_20181116_0037353.mat'
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
