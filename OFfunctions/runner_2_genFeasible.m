%%%find nearest feasible soln based on random/user parameters%%%
%%%choose type of simulation: SBR-mode or ODE-mode %%%
%%%no data required, to test whether model works%%%
%%%simulation of active pool only, backbone carbon only%%%
%%%exclude stagnant pool (which is data dependent)%%%

OFspec = OpenFLUX.OFobjSpecification(mfilename,OFspecFileName);

load(OFspec.modelObjSaveName);

if ~isfolder(OFspec.simSaveFolder)
    mkdir(OFspec.simSaveFolder);
end

OF.genLabelledSubstrate;
simParas = OF.prepSimulation;

opOptions = optimoptions('fmincon', 'Display','iter','MaxFunEvals',40000);

x0 = simParas.lb + rand(size(simParas.lb)).*(simParas.ub-simParas.lb);
xGuess = x0;

%%%find a feasible solution
if OF.isODEsolver%%ODE-mode
    while 1
        xFeas = fmincon(@(x)minDistX0(x,x0),xGuess,[],[],[],[],simParas.lb,simParas.ub,simParas.conFxn,opOptions);
        if all(simParas.conFxn(xFeas)<=0)
            break
        else
            xGuess = xFeas;
        end
    end
    if OFspec.bumpUpXFeasIniConc
        xFeas = OF.simSolnODE_stepConc(xFeas);
    end    
elseif OF.isDynamic%%%SBR-mode
    while 1
        xFeas = fmincon(@(x)minDistX0(x,x0),xGuess,simParas.Acon,simParas.Bcon,[],[],simParas.lb,simParas.ub,[],opOptions);
        if all(simParas.Acon*xFeas<=simParas.Bcon)
            break
        else
            xGuess = xFeas;
        end
    end
else%%SS mode
    while 1
        xFeas = fmincon(@(x)minDistX0(x,x0),xGuess,[],[],[],[],simParas.lb,simParas.ub,simParas.conFxn,opOptions);
        if all(simParas.conFxn(xFeas)<=0)
            break
        else
            xGuess = xFeas;
        end
    end
    
    %min-max simulate
    opOptions = optimoptions('fmincon', 'Display','off','MaxFunEvals',40000);
    fluxRange = zeros(size(OF.opInput.ns_free,1),2);
    for i = 1:size(OF.opInput.ns_free,1)
        vChoose = -i;
        fOBJ = @(x)fluxSimRange(x,OF.opInput.pEntries,OF.opInput.isRevPivot,OF.opInput.par,OF.opInput.ns_free,OF.opInput.v_fixed,vChoose);
        xMax = fmincon(fOBJ,xFeas,[],[],[],[],simParas.lb,simParas.ub,simParas.conFxn,opOptions);
        vChoose = i;
        fOBJ = @(x)fluxSimRange(x,OF.opInput.pEntries,OF.opInput.isRevPivot,OF.opInput.par,OF.opInput.ns_free,OF.opInput.v_fixed,vChoose);
        xMin = fmincon(fOBJ,xFeas,[],[],[],[],simParas.lb,simParas.ub,simParas.conFxn,opOptions);
        fluxRange(i,:) = [fOBJ(xMin) fOBJ(xMax)];
    end
end


simSave.x0 = x0;
simSave.xFeas = xFeas;
simSave.fluxRange = fluxRange;
simSave.datetimeCreated = datetime('now','Format','yyyyMMdd_HHmmSSS');
if OF.isODEsolver%%ODE-mode
    simSave.saveFileName = strcat(['ODEsim_' char(simSave.datetimeCreated) '.mat']);
elseif OF.isDynamic%%%SBR-mode
    simSave.saveFileName = strcat(['SBRsim_' char(simSave.datetimeCreated) '.mat']);
else %%%SS-mode
    simSave.saveFileName = strcat(['SSsim_' char(simSave.datetimeCreated) '.mat']);
end
save(strcat(OFspec.simSaveFolder, simSave.saveFileName), 'simSave','OF');
disp(strcat(OFspec.simSaveFolder, simSave.saveFileName));

if OFspec.simulateSoln 
    simOutput = OF.simSoln(xFeas);
    save(strcat(OFspec.simSaveFolder, simSave.saveFileName), 'simSave','OF','simOutput');
    disp('simOutput generated');
end

clearvars -except OF OFspec simOutput simSave