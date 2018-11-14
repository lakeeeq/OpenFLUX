%%%run simulation based on random/user parameters%%%
%%%no data required, to test whether model works%%%
%%%simulation of active pool only, backbone carbon only%%%
%%%exclude stagnant (data dependent)%%%
clearvars
clc
reLoadData = true;
% reLoadData = false;
if reLoadData
    load OFobj_ODEsim
else
    load OFobj
    OF.isODEsolver = true;
    OF.genLabelledSubstrate;
    simParas = OF.prepSimulation;
    save OFobj_ODEsim OF simParas
end
x0 = rand(size(simParas.lb));
load xFeasSBRsim
% x0 = xFeas;
% opOptions = optimoptions('fmincon', 'Display','iter','MaxFunEvals',40000);
% 
% simParas.conFxn(x0);
% % return
% %%%find a feasible solution
% xGuess = x0;
% while 1
%     xFeas = fmincon(@(x)minDistX0(x,x0),xGuess,[],[],[],[],simParas.lb,simParas.ub,simParas.conFxn,opOptions);
%     if all(simParas.conFxn(xFeas)<=0)
%         break
%     else
%         xGuess = xFeas;
%     end
% end
return
[simEMU,simConc,simFlux,simTime] = OF.simSoln(xFeas);

