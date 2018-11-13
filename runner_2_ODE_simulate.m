%%%run simulation based on random/user parameters%%%
%%%no data required, to test whether model works%%%%
%%%simulation of active pool only, exclude stagnant (data dependent)%%%%
clearvars
clc
load OFobj
[opCon,opInput] = OF.prepSimulation;
x0 = rand(size(opCon.lb));
opOptions = optimoptions('fmincon', 'Display','iter','MaxFunEvals',40000);

%%%find a feasible solution
xGuess = x0;
while 1
    xFeas = fmincon(@(x)minDistX0(x,x0),xGuess,opCon.Acon,opCon.Bcon,[],[],opCon.lb,opCon.ub,[],opOptions);
    if all(opCon.Acon*xFeas<=opCon.Bcon)
        break
    else
        xGuess = xFeas;
    end
end

[simEMU,simConc,simFlux,simTime] = simSoln(OF,opInput,xFeas);

