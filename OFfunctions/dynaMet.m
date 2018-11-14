function dy = dynaMet(t,y,y_empty,knotSeq,orderS,noEMUperm,...
    odeModel,EMUstate_struct,cauchyTags,fluxStoicT,CP,Y2concY,maxT)
%removed Y2conc, conc2X
dy = y_empty;
x = y./(Y2concY*y);

%%%%repopulate EMUstate matrix
NoutCol = bSplineMat_lite(knotSeq,t/maxT,orderS);
fluxCombine = fluxStoicT*CP*NoutCol;%%%flux needs int

for i = 1:noEMUperm
    %%% step 1%% do calc matrix
    if odeModel(i).calcSize == 1
        EMUstate_struct(i).calcEMU = x(odeModel(i).y2isoRows_Shaped)';
    else
        EMUstate_struct(i).calcEMU = x(odeModel(i).y2isoRows_Shaped);
    end
    
    %%%don't need to change input substrate rows
    matNow = cauchyTags{i};
    for k = 1:size(matNow,1)
        state1 = EMUstate_struct(matNow(k,1)).calcEMU(matNow(k,2),:);
        state2 = EMUstate_struct(matNow(k,3)).calcEMU(matNow(k,4),:);
        EMUstate_struct(i).cauchyEMU(k,:) = cauchy(state1,state2);
    end
    A = odeModel(i).A0;
    A(odeModel(i).Amap(:,1)) = fluxCombine(odeModel(i).Amap(:,2));
    dy_mat = A * [EMUstate_struct(i).calcEMU; EMUstate_struct(i).isEMU; EMUstate_struct(i).cauchyEMU];
    dy(odeModel(i).y2isoRows) = dy_mat;
end
% dy(y<0) = 1e-10;