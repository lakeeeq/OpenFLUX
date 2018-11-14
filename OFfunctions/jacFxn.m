function Jmat = jacFxn(t,Y_in,jacOut,knotSeq,CP,maxT)
NoutCol = bSplineMat_lite(knotSeq,t/maxT,jacOut.orderS);
fluxCombine = jacOut.fluxStoicTode*CP*NoutCol;%%%flux needs int

k = 1./(jacOut.Y2concY*Y_in);
xFract = Y_in.*k;
xCoeff1 = k - Y_in.*k.^2;
xCoeff2 = -xFract(jacOut.yTagBEM_denom).*k(jacOut.yTagBEM_denom);%col2 coeff on top, col 1 denominator
cc1 = k(jacOut.cauchyCoeff1);
cc2 = jacOut.cc2;
cc2(jacOut.nzEle_cc2) = xFract(jacOut.cauchyCoeff2(jacOut.nzEle_cc2));
cc3DP1 = jacOut.cc3DP1;
cc3DP2 = jacOut.cc3DP2;
cc3DP1(jacOut.nzEle_cc3) = xFract(jacOut.cc3DP1(jacOut.nzEle_cc3));
cc3DP2(jacOut.nzEle_cc3) = xFract(jacOut.cc3DP2(jacOut.nzEle_cc3));
cc3 = sum(cc3DP1.*cc3DP2,2);
xCauchy = cc1.*(cc2-cc3);
xCoeffCombine = [xCoeff1;xCoeff2;xCauchy];

Jcoeff_0 = jacOut.Jcoeff_0;
Jcoeff_0(jacOut.Jcoeff_map(:,1)) = xCoeffCombine(jacOut.Jcoeff_map(:,2));
Jvect = Jcoeff_0*fluxCombine;
Jmat = jacOut.Jmat;
Jmat(jacOut.JnzRows) = Jvect;