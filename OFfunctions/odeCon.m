function [c,ceq] = odeCon(x,c_base,concMap,c_diff,f_base,CPmap,f_diff,...
    Sint,knotSeq,orderS,tSimStep,maxT)

Y0 = c_base + x(concMap).*c_diff;
CP = f_base + x(CPmap).*f_diff;

odeFxn = @(t,y)odeMet(t,y);

% [t,y] = ode45(odeFxn,tSimStep,Y0);
[t,y] = ode15s(odeFxn,tSimStep,Y0);
ceq = 0;
c = -y(:)+0.01;

    function dy = odeMet(t,y)
        NoutCol = bSplineMat_lite(knotSeq,t/maxT,orderS);
        dy = Sint*CP*NoutCol;
    end
end