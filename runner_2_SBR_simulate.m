%%%run simulation based on random/user parameters%%%
%%%no data required, to test whether model works%%%
%%%simulation of active pool only, backbone carbon only%%%
%%%exclude stagnant (data dependent)%%%
%%%knot position fixed to settings%%%
clearvars
clc
% rng('shuffle');
reLoadData = true;
reLoadData = false;
if reLoadData
    load SBR_OFobj_sim
else
    load OFobj
    OF.isODEsolver = false;
    OF.genLabelledSubstrate;
    simParas = OF.prepSimulation;
    save SBR_OFobj_sim OF simParas
end

x0 = rand(size(simParas.lb));
opOptions = optimoptions('fmincon', 'Display','iter','MaxFunEvals',40000);

%%%find a feasible solution
xGuess = x0;
while 1
    xFeas = fmincon(@(x)minDistX0(x,x0),xGuess,simParas.Acon,simParas.Bcon,[],[],simParas.lb,simParas.ub,[],opOptions);
    if all(simParas.Acon*xFeas<=simParas.Bcon)
        break
    else
        xGuess = xFeas;
    end
end

[simEMU,simConc,simFlux,simTime] = OF.simSoln(xFeas);