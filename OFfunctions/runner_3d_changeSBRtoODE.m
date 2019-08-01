%%%%change SBR to ODE, and vice versa%%%%
%%%create from scratch, but copy over knot position and existing solution
%%%load and change specific (listed) optimisation instances%%%

OFspec = OpenFLUX.OFobjSpecification(mfilename,OFspecFileName);
load(OFspec.modelObjSaveName);
OF.isODEsolver = OFspec.isODEsolver;%change from SBR to ODE solver, and vice versa
OF.isOptimisation = true;
OF.genLabelledSubstrate;

%specify extra constraints
OFadditionalConstraintsScript

if OF.isODEsolver %SBR to ODE
    OF.odeSimTime = OFspec.odeSimTime;
    for i = 1:numel(OFspec.fileListToChange)
        disp(['loading ' OFspec.fileListToChange{i}]);
        load(strcat([OFspec.opSaveFolder OFspec.fileListToChange{i}]));
        if opSave.isODE
            disp('this instance is already an ODE. Skip')
            continue
        end
        OF.metDataFileName = opSave.metDataFileName;
        OF.reEstimateError('load',[]);
        OF.intKntPos = opSave.xIntKnot;
        simParas = OF.prepSimulation;
        if i == 1
            disp('jacobian is generated for the first instance, then copied to the rest')
        end
        
        %%%%%%%%%%%
        %specify extra data (radiation)
        OFadditionalDataScript
        %%%%%%%%%%%
        
        fitFxn = OF.generateFitFxn(additionalData);
        opSave.fitFxn = fitFxn;
        opSave.conFxn = simParas.conFxn;
        opSave = rmfield(opSave,'AconParas');
        opSave = rmfield(opSave,'stepBTWsample');
        opSave.isODE = true;
        opSave.saveFileName = strcat(['ODEop_' char(opSave.datetimeCreated) '.mat']);
        
        save(strcat(OFspec.opSaveFolder, opSave.saveFileName), 'opSave');
        disp(['created ' opSave.saveFileName]);
        fprintf('\n');
    end
    
    
else %ODE to SBR
    OF.odeSimTime = [];
    for i = 1:numel(OFspec.fileListToChange)
        disp(['loading ' OFspec.fileListToChange{i}]);
        load(strcat([OFspec.opSaveFolder OFspec.fileListToChange{i}]));
        if ~opSave.isODE
            disp('this instance is already an SBR. Skip')
            continue
        end
        OF.metDataFileName = opSave.metDataFileName;
        OF.reEstimateError('load',[]);
        OF.intKntPos = opSave.xIntKnot;
        simParas = OF.prepSimulation;
       
        
        %%%%%%%%%%%
        %specify extra data (radiation)
        OFadditionalDataScript
        %%%%%%%%%%%
        
        opSave.stepBTWsample = OFspec.stepBTWsample;
        opSave.AconParas = simParas.AconParas;
        
        fitFxn = OF.generateFitFxn(additionalData);
        opSave.fitFxn = fitFxn;
        opSave = rmfield(opSave,'conFxn');
        opSave.isODE = false;
        opSave.saveFileName = strcat(['SBRop_' char(opSave.datetimeCreated) '.mat']);
        
        save(strcat(OFspec.opSaveFolder, opSave.saveFileName), 'opSave');
        disp(['created ' opSave.saveFileName]);
        fprintf('\n');
    end    
end