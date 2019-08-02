%%%%change SBR to ODE, and vice versa%%%%
%%%create from scratch, but copy over knot position and existing solution
%%%load and change specific (listed) optimisation instances%%%

OFspec = OpenFLUX.OFobjSpecification(mfilename,OFspecFileName);
% load(OFspec.modelObjSaveName);
% OF.isODEsolver = OFspec.isODEsolver;%change from SBR to ODE solver, and vice versa
% OF.isOptimisation = true;
% OF.genLabelledSubstrate;

%specify extra constraints
% OFadditionalConstraintsScript

% if OF.isODEsolver %change SBR to ODE
%     OF.odeSimTime = OFspec.odeSimTime;

if OFspec.isODEsolver
    disp('changing from SBR to ODE');
else
    disp('changing from ODE to SBR');
end
doneOnce = false;
for i = 1:numel(OFspec.fileListToChange)
    disp(['loading ' OFspec.fileListToChange{i}]);
    try
        load(strcat([OFspec.opSaveFolder OFspec.fileListToChange{i}]));
    catch
        disp(['could not find/load ' OFspec.fileListToChange{i}]);
        continue
    end
    if opSave.isODE && OFspec.isODEsolver
        disp('this instance is already an ODE. Skip')
        continue
    end
    if ~opSave.isODE && ~OFspec.isODEsolver
        disp('this instance is already an SBR. Skip')
        continue
    end
    OF.isODEsolver = OFspec.isODEsolver;
    
    if OFspec.isODEsolver%%%SBR to ODE
        OF.odeSimTime = OFspec.odeSimTime;
        if doneOnce
          OF.writePrivProp('jacOut',newJacOut);  
        end
        simParas = OF.prepSimulation;
        if ~doneOnce
            disp('jacobian is generated for the first instance, then copied to the rest')
            newJacOut = OF.getPrivProp('jacOut');
            doneOnce = true;
        end
        fitFxn = OF.generateFitFxn;
        opSave.fitFxn = fitFxn;
        opSave.conFxn = simParas.conFxn;
        opSave = rmfield(opSave,'AconParas');
        opSave = rmfield(opSave,'stepBTWsample');
        opSave.isODE = true;
        opSave.saveFileName = strcat(['ODEop_' char(opSave.datetimeCreated) '.mat']);
                
    else%%%ODE to SBR
        OF.odeSimTime = [];
        OF.stepBTWsample = OFspec.newSampleSteps;
        simParas = OF.prepSimulation;        
        opSave.AconParas = simParas.AconParas;
        
        fitFxn = OF.generateFitFxn;
        opSave.fitFxn = fitFxn;
        opSave.stepBTWsample = OF.stepBTWsample;
        opSave = rmfield(opSave,'conFxn');
        opSave.isODE = false;
        opSave.saveFileName = strcat(['SBRop_' char(opSave.datetimeCreated) '.mat']);

    end
    save(strcat(OFspec.opSaveFolder, opSave.saveFileName), 'opSave', 'OF');
    disp(['created ' opSave.saveFileName]);
    fprintf('\n');
end
