function [c,ceq] = odeCon(x,c_base,concMap,c_diff,f_base,CPmap,f_diff,...
    Sint,knotSeq,orderS,tSimStep,maxT,minConc)

Y0 = c_base + x(concMap).*c_diff;
CP = f_base + x(CPmap).*f_diff;

odeFxn = @(t,y)odeMet(t,y,knotSeq,maxT,orderS,Sint,CP);

[t,y] = ode15s(odeFxn,tSimStep,Y0);
ceq = 0;
c = -y(:)+minConc;
