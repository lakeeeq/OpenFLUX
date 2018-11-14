%%%create OpenFLUX object and read txt model%%%
%%%specify input substrate%%%
%%%specify time of samples%%%
%%%specifcy ion formula file%%% 
clearvars
clc
addpath ./OFfunctions;
addpath ./inputs
OF = OpenFLUX;
% return
%%%%initial configuration%%%%
OF.modelCondition = 'insulin';%tag for OF object
OF.modelFileName = 'model.txt';%metabolic model
OF.ionFormFileName = 'metIonFormula.txt';%formula for MZs
OF.buildModel;

OF.intKntPos = [0.2];%specify bspline internal knots
OF.orderS = 3;%specify bspline order

OF.labelledSub = {'GLC_in' [0.99 0.99 0.99 0.99 0.99 0.99]};%%%specify positional enrichment
OF.sampleTime = [1 5 10 20 40 60];%indicate sampling time points (still required for simulations)

% OF.stepBTWsample = [20	20	20	20	20	20];%required for SBR
OF.stepBTWsample = [500	400	200	75	75	50];%required for SBR
OF.odeSimTime = [0:1:60];%required for ODE15s, sampleTime must be within this range

save OFobj OF
