function [f, metConcProfile,EMUstateStore,vProfile,fConc_diff,fMID_diff,rad_diff,vInt] =...
    leastSQ_SBR(x,CPmap,concMap,f_base,f_diff,c_base,c_diff,noT,Nout_int,...
    noEMUperm,cauchyTags,Sint,deltaT,tSampleIndex,noTsample,dataMet,fluxStoicT,...
    Nout,noExpData,dataMetConc_EXPvect,dataMetSE_EXPvect,dataMetMID_SIMvect,...
    dataMetMID_EXPvect,dataMetMID_SEvect,EMUstateStoreIS_block,EMUstateStore,...
    EMUstate,stagMap,cstag,metConcProfile_SIM,A_cell,Cmap_cell,Vmap_cell,EMUsize,...
    additionalData)

%%Forward Euler %EMUstate_i+1 = A*EMUstate_i
CP = f_base + x(CPmap).*f_diff;
vProfile = CP*Nout;
metConc = c_base + x(concMap).*c_diff;%this is the active pool
stagAmount = x(stagMap).*cstag;%active fraction, range is 0.001 to 1
vInt = CP*Nout_int;
metConcIni = metConc(:,ones(noT,1));
metConcProfile = Sint*vInt + metConcIni;

EMUstateStore{1} = EMUstate(:,1);
for i = 2:noT
    flux = vProfile(:,i)*deltaT(i);
    fluxCombine = fluxStoicT*flux;%convert individual fluxes to summary parameters
    
    %update EMU state
    for j = noEMUperm:-1:1
        A = A_cell{j};%empty A matrix
        Cmap = Cmap_cell{j};
        Vmap = Vmap_cell{j};
        A(Cmap(:,1)) = metConc(Cmap(:,2));
        A(Vmap(:,1)) = fluxCombine(Vmap(:,2));
        EMUstateNext = A*[EMUstate{j,1};EMUstate{j,2}];%individual top, cauchy bottom
        EMUstateNextRowSum = sum(EMUstateNext,2);
        EMUstateNext = EMUstateNext./EMUstateNextRowSum(:,ones(1,EMUsize(j)));
        EMUstate{j,1} = [EMUstateNext; EMUstateStoreIS_block{i,j}];%update individual only
    end
    
    %update EMU cauchy state, for next round
    for j = 1:noEMUperm
        matNow = cauchyTags{j};
        for k = 1:size(matNow,1)
            state1 = EMUstate{matNow(k,1),1}(matNow(k,2),:);
            state2 = EMUstate{matNow(k,3),1}(matNow(k,4),:);
            cauchyProd = cauchy(state1,state2);
            EMUstate{j,2}(k,:) = cauchyProd;
        end
    end
    EMUstateStore{i} = EMUstate(:,1);
    metConc = metConcProfile(:,i);
end

noTSample = 1:noTsample;
for i = 1:noExpData
    midSize = dataMet{i,12};
    mid_SIM = dataMet{i,6};%clone
    hitX = dataMet{i,10};
    for j = noTSample
        mid_SIM(:,j) = EMUstateStore{tSampleIndex(j)}{hitX(1)}(hitX(2),:);
    end
    
    %%add stagnant unlabeled pool
    activAmtMat = metConcProfile(dataMet{i,8}*ones(midSize,1),tSampleIndex);
    midSIM_new = mid_SIM.*activAmtMat + dataMet{i,13}*stagAmount(i);    

    midCorr_SIM = dataMet{i,11}*midSIM_new;
    midCorr_SIM(dataMet{i,14}) = 0;
    midCorr_SIM_sum = sum(midCorr_SIM,1);
    midCorr_SIM = midCorr_SIM./midCorr_SIM_sum(ones(midSize,1),:);
    dataMetMID_SIMvect(dataMet{i,9}(:,2)) = midCorr_SIM(dataMet{i,9}(:,1));
    metConcProfile_SIM(i,:) = midCorr_SIM_sum;
end
metConcProfile_SIMvect = metConcProfile_SIM(:);


fConc_diff = (dataMetConc_EXPvect-metConcProfile_SIMvect)./dataMetSE_EXPvect;
fMID_diff = (dataMetMID_EXPvect-dataMetMID_SIMvect)./dataMetMID_SEvect;
f = fConc_diff'*fConc_diff + fMID_diff'*fMID_diff;
