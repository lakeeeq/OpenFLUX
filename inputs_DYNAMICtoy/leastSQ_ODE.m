function [f, Yj,fConc_diff,fMID_diff,rad_diff] = leastSQ_ODE(x,c_base,concMap,c_diff,conc2X,Y0f,f_base,CPmap,...
    f_diff,yLength,knotSeq,orderS,noEMUperm,odeModel,...
    EMUstate_struct,cauchyTags,fluxStoicTode,Y2concY,jacOut,tSample,...
    mFmap,mCmap,Ymat_Czeroed,mFmap_concDenon,dataMetConc_EXPvect,dataMetMID_EXPvect,...
    dataMetSE_EXPvect,dataMetMID_SEvect,stagMap,cstag,stag_2_Y_map,Y_CM_stag,...
    YstagVect,YsimVect,dataMet,additionalData,noExpData,maxT)

metConcIni = c_base + x(concMap).*c_diff;
Y0 = metConcIni(conc2X).*Y0f;
stagAmount = x(stagMap).*cstag;%active fraction, range is 0.001 to 1

CP = f_base + x(CPmap).*f_diff;

y_empty = zeros(yLength,1);
dynaFxn = @(t,y)dynaMet(t,y,y_empty,knotSeq,orderS,noEMUperm,...
    odeModel,EMUstate_struct,cauchyTags,fluxStoicTode,CP,Y2concY,maxT);

dynaJacFxn = @(t,y)jacFxn(t,y,jacOut,knotSeq,CP,maxT);
odeOptionsJ = odeset('NonNegative',1,'Jacobian',dynaJacFxn,'RelTol',1e-3);
while 1
    [Tj, Yj] = ode15s(dynaFxn,[0 tSample],Y0,odeOptionsJ);
    if numel(Tj) < numel(tSample)+1
        odeOptionsJ.RelTol =  odeOptionsJ.RelTol*10;
        disp(['err, repeat with larger 10x reltol to ' num2str(odeOptionsJ.RelTol)]);
    else
        break
    end
end

%%%%factoring in stagnant pool
Ycopy = Yj(2:end,:)';
YstagVect(stag_2_Y_map(:,1)) = stagAmount(stag_2_Y_map(:,2));
Ycopy2 = Ycopy+Y_CM_stag.*YstagVect;

%%%%correcting for natural enrichment
for i = 1:noExpData
    Ycopy2(dataMet{i,16},:) = dataMet{i,11}*Ycopy2(dataMet{i,16},:);
end

%%%make sure concentrations are summed from the corrected MIDs
Ycopy2(Ymat_Czeroed) = 0;
YsimVect(mFmap(:,2)) = Ycopy2(mFmap(:,1));
CsimVect = mCmap*Ycopy2;
CsimVect_15sJ = CsimVect(:);
YsimVect_15sJ = YsimVect./CsimVect(mFmap_concDenon);


fConc_diff = (dataMetConc_EXPvect-CsimVect_15sJ)./dataMetSE_EXPvect;
fMID_diff = (dataMetMID_EXPvect-YsimVect_15sJ)./dataMetMID_SEvect;

f = fConc_diff'*fConc_diff + fMID_diff'*fMID_diff;
rad_diff=[];