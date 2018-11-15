%%%run simulation based on random/user parameters%%%
%%%no data required, to test whether model works%%%
%%%simulation of active pool only, backbone carbon only%%%
%%%exclude stagnant (data dependent)%%%
clearvars
clc
% rng('shuffle');
reLoadData = true;
reLoadData = false;
if reLoadData
    load ODE_OFobj_sim
else
    load OFobj
    OF.isODEsolver = true;
    OF.genLabelledSubstrate;
    simParas = OF.prepSimulation;
    save ODE_OFobj_sim OF simParas
end
x0 = rand(size(simParas.lb));
opOptions = optimoptions('fmincon', 'Display','iter','MaxFunEvals',40000);

simParas.conFxn(x0);

%%%find a feasible solution
xGuess = x0;
while 1
    xFeas = fmincon(@(x)minDistX0(x,x0),xGuess,[],[],[],[],simParas.lb,simParas.ub,simParas.conFxn,opOptions);
    if all(simParas.conFxn(xFeas)<=0)
        break
    else
        xGuess = xFeas;
    end
end

[simEMU,simConc,simFlux,simTime] = OF.simSoln(xFeas);