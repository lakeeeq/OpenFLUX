function [f, EMUstate,v,fMID_diff] =...
    leastSQ_SS(x,EMUstate,pEntries,isRevPivot,par,ns_free,v_fixed,fluxStoicT,...
    A_cell,B_cell,Amap,Bmap,cauchyMap,knownMapLHS,calcMapLHS,knownMap,calcMap,noEMUperm,...
    noExpData,dataMet,dataMetMID_SIMvect,dataMetMID_EXPvect,dataMetMID_SEvect,...
    additionalData)


xFlux = x(~pEntries);
pFract = x(pEntries);
xFlux(isRevPivot) = my_hyperfunc(xFlux(isRevPivot),par);
v = ns_free*xFlux + v_fixed;
fluxCombine = fluxStoicT*v;
for i = 1:noEMUperm
    A = A_cell{i};
    B = B_cell{i};
    A(Amap{i}(:,1)) = fluxCombine(Amap{i}(:,2));
    B(Bmap{i}(:,1)) = fluxCombine(Bmap{i}(:,2));
    %do cauchy list
    cauchyInput = [];
    for j = 1:size(cauchyMap{i},1)
        state1 = EMUstate{cauchyMap{i}(j,1)}(cauchyMap{i}(j,2),:);
        state2 = EMUstate{cauchyMap{i}(j,3)}(cauchyMap{i}(j,4),:);
        cauchyInput(j,:) = cauchy(state1,state2);
    end
    
    if knownMapLHS(i)==0
        Y = cauchyInput;
    else
        Y = [EMUstate{knownMapLHS(i)}(knownMap{i},:); cauchyInput];
    end
    
    X = A\B*Y;
    EMUstate{calcMapLHS(i)}(calcMap{i},:) = X;
end

for i = 1:noExpData
    midSize = dataMet{i,12};
    hitX = dataMet{i,10};
    mid_SIM = EMUstate{hitX(1)}(hitX(2),:);
    
    %%add ISA unlabeled pool
    if dataMet{i,15}~=0
        mid_SIM = mid_SIM*(1-pFract(dataMet{i,15})) + dataMet{i,13}*pFract(dataMet{i,15});
    end
    
    midCorr_SIM = dataMet{i,11}*mid_SIM';
    midCorr_SIM(dataMet{i,14}) = 0;
    %     midCorr_SIM_sum = sum(midCorr_SIM,1);
    midCorr_SIM = midCorr_SIM/sum(midCorr_SIM);
    dataMetMID_SIMvect(dataMet{i,9}(:,2)) = midCorr_SIM(dataMet{i,9}(:,1));
end

fMID_diff = (dataMetMID_EXPvect-dataMetMID_SIMvect)./dataMetMID_SEvect;
f = fMID_diff'*fMID_diff;