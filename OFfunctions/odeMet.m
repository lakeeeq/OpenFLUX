function dy = odeMet(t,y,knotSeq,maxT,orderS,Sint,CP)
NoutCol = bSplineMat_lite(knotSeq,t/maxT,orderS);
dy = Sint*CP*NoutCol;
