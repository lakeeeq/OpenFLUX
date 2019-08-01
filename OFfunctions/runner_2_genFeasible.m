%%%find nearest feasible soln based on random/user parameters%%%
%%%choose type of simulation: SBR-mode or ODE-mode %%%
%%%no data required, to test whether model works%%%
%%%simulation of active pool only, backbone carbon only%%%
%%%exclude stagnant pool (which is data dependent)%%%

OFspec = OpenFLUX.OFobjSpecification(mfilename,OFspecFileName);

load(OFspec.modelObjSaveName);
OF.genLabelledSubstrate;
simParas = OF.prepSimulation;

opOptions = optimoptions('fmincon', 'Display','iter','MaxFunEvals',40000);

x0 = rand(size(simParas.lb));
xGuess = x0;

%%%find a feasible solution
if OF.isODEsolver%%ODE-mode
    simParas.conFxn(x0);
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
    
else%%%SBR-mode
    while 1
        xFeas = fmincon(@(x)minDistX0(x,x0),xGuess,simParas.Acon,simParas.Bcon,[],[],simParas.lb,simParas.ub,[],opOptions);
        if all(simParas.Acon*xFeas<=simParas.Bcon)
            break
        else
            xGuess = xFeas;
        end
    end
end

simSave.x0 = x0;
simSave.xFeas = xFeas;
simSave.datetimeCreated = datetime('now','Format','yyyyMMdd_HHmmSSS');
if OF.isODEsolver%%ODE-mode
    simSave.saveFileName = strcat(['ODEsim_' char(simSave.datetimeCreated) '.mat']);
else
    simSave.saveFileName = strcat(['SBRsim_' char(simSave.datetimeCreated) '.mat']);
end
save(strcat(OFspec.simSaveFolder, simSave.saveFileName), 'simSave','OF');
disp(strcat(OFspec.simSaveFolder, simSave.saveFileName));

% [simEMU,simConc,simFlux,simTime] = OF.simSoln(xFeas);