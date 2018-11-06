%%%create OpenFLUX object and read txt model%%%
%%%specify input substrate%%%
%%%specify time of samples%%%
%%%specifcy ion formula file%%% 
clearvars
clc
addpath ./OFfxns;
OF = OpenFLUX;
% return
%%%%initial configuration%%%%
OF.modelCondition = 'insulin';
OF.modelFileName = 'model.txt';
OF.buildModel;
labSub = {'GLC_in' [0.99 0.99 0.99 0.99 0.99 0.99]};%%%specify positional enrichment
OF.genLabelledSubstrate(labSub);
OF.sampleTime = [1 5 10 20 40 60];
% OF.stepBTWsample = [20	20	20	20	20	20];
OF.stepBTWsample = [500	400	200	75	75	50];
OF.intKntPos = 0.2;
OF.ionFormFileName = 'metIonFormula.txt';

save OFobj OF
