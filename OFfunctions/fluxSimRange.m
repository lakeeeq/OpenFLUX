function f = fluxSimRange(x,pEntries,isRevPivot,par,ns_free,v_fixed,vChoose)
%for steady state
xFlux = x(~pEntries);
xFlux(isRevPivot) = my_hyperfunc(xFlux(isRevPivot),par);
v = ns_free*xFlux + v_fixed;
f = sign(vChoose)*v(abs(vChoose));
