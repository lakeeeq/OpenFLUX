function [cIneq,cEq] = conFluxSim(x,pEntries,isRevPivot,par,ns_free,v_fixed)

xFlux = x(~pEntries);
xFlux(isRevPivot) = my_hyperfunc(xFlux(isRevPivot),par);
v = ns_free*xFlux + v_fixed;
cIneq = -v+0.001;
cEq=[];