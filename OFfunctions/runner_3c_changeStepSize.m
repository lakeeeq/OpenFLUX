%%%%change step size of SBR%%%
%%%load and change specific (listed) optimisation instances%%%
% clearvars
% clc
OFspec = OpenFLUX.OFobjSpecification(mfilename,OFspecFileName);

% load(OFspec.modelObjSaveName);
% if OF.isODEsolver
%     disp('step size change only applies to SBR')
%     return
% end
% OF.isOptimisation = true;
% OF.genLabelledSubstrate;
% 
% %specify extra constraints
% OFadditionalConstraintsScript
disp(['changing step size to: ' num2str(OFspec.newSampleSteps)]);
for i = 1:numel(OFspec.fileListToChange)
    try
        load(strcat([OFspec.opSaveFolder OFspec.fileListToChange{i}]));
    catch
        disp(['could not find/load ' OFspec.fileListToChange{i}]);
        continue
    end
    if opSave.isODE
        disp(['this op instance is ODE, skip: ' OFspec.fileListToChange{i}]);
        continue
    else
        disp(['loading ' OFspec.fileListToChange{i}]);
    end
%     OF.metDataFileName = opSave.metDataFileName;%for insulin
%     OF.reEstimateError('load',[]);
    disp(['original step size: ' num2str(opSave.stepBTWsample)]);
    OF.stepBTWsample = OFspec.newSampleSteps;
%     OF.intKntPos = opSave.xIntKnot;
    simParas = OF.prepSimulation;
       
    %%%%%%%%%%%
    %specify extra data (radiation)
%     OFadditionalDataScript
    %%%%%%%%%%%
    
    fitFxn = OF.generateFitFxn;
    opSave.fitFxn = fitFxn;
    opSave.AconParas = simParas.AconParas;
    opSave.stepBTWsample = OF.stepBTWsample;
    save(strcat(OFspec.opSaveFolder, opSave.saveFileName), 'opSave', 'OF');
end
