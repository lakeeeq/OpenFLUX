clearvars
clc
addpath ./OFfxns;
OF = OpenFLUX;
% return
%%%%initial configuration%%%%
OF.modelCondition = 'adipocyte insulin 1 hr';
OF.modelFileName = 'model.txt';
OF.buildModel;
labSub = {'GLC_in' [0.99 0.99 0.99 0.99 0.99 0.99]
    'GLNx_in' [0.99 0.0107 0.0107 0.0107 0.0107]};%%%specify positional enrichment
OF.genLabelledSubstrate(labSub);
OF.sampleTime = [1 5 10 20 40 60];
OF.stepBTWsample = [20	20	20	20	20	20];
% OF.stepBWTsample = [500	400	200	75	75	50];
OF.intKntPos = 0.2;
OF.ionFormFileName = 'metIonFormula.txt';

%%%%setup data input%%%%
OF.metDataFileName = 'metDataIns.txt';
% OF.importMetData;
% OF.reEstimateError('generate',10000);
% OF.reEstimateError('save',[],'metData_insulin');
OF.reEstimateError('load',[],'metData_insulin');

%%%%manual modifications%%%%
%constrain accoa_out and glycogen_out initial to near zero 0.01
hitMet = strcmp('ACCOA_out',OF.metListInt);
OF.concScale(hitMet,:) = 0.01;
hitMet = strcmp('GLYCOGEN_out',OF.metListInt);
OF.concScale(hitMet,:) = 0.01;


%%%%setup simulation%%%%
[opCon,opInput] = OF.prepSimulation;
x0 = rand(size(opCon.lb));
opOptions = optimoptions('fmincon', 'Display','iter','MaxFunEvals',40000);

%{
xGuess = x0;
while 1
    xFeas = fmincon(@(x)minDistX0(x,x0),xGuess,opCon.Acon,opCon.Bcon,[],[],opCon.lb,opCon.ub,[],opOptions);
    if all(opCon.Acon*xFeas<=opCon.Bcon)
        break
    else
        xGuess = xFeas;
    end
end
save xFeas xFeas
return
%}
load xFeas

%%radiolabel data
%%% FA accumulation data
%c1 glucose, c2 mix, c3 unlabeled
radData{1,1} = [cVectGen2([0.01 0.99],2,3)'
    cauchy(cVectGen2([1-0.0107 0.0107],1,2)',cVectGen2([0.01 0.99],1,2)')
cVectGen2([1-0.0107 0.0107],2,3)']';
radData{1,1} = pinv(radData{1,1});

radData{2,1} = [cVectGen2([0.01 0.99],6,7) cVectGen2([0.9893 0.0107],6,7)];
radData{2,1} = pinv(radData{2,1});

%%FA accumulation pCmol/mgP accumulated for 1 hr
if strcmp('metDataIns.txt',OF.metDataFileName)
    radData{1,2} = [76973.028 9680.29035];%average ,se for FA
    radData{2,2} = [241.2611 10.16552]*6000;%average ,se for glycogen
else
   radData{1,2} = [16819.87	652.2888155];%average ,se for FA
   radData{2,2} = [6.413856	0.478499]*6000;%average ,se for glycogen
end
bigEMUmodel = getPrivProp(OF,'bigEMUmodel');
hitRow = cell2mat(bigEMUmodel(:,1))==2;
hitEMU = matchEMU('ACCOA_out',[1 1],bigEMUmodel{hitRow,2});
hitIntMet = strcmp('ACCOA_out',OF.metListInt);
radData{1,3} = [find(hitRow) find(hitEMU) find(hitIntMet)];
hitRow = cell2mat(bigEMUmodel(:,1))==6;
hitEMU = matchEMU('GLYCOGEN_out',[1 1 1 1 1 1],bigEMUmodel{hitRow,2});
hitIntMet = strcmp('GLYCOGEN_out',OF.metListInt);
radData{2,3} = [find(hitRow) find(hitEMU) find(hitIntMet)];
%%%%%%%%%%%%%%%%%%%%%%


%%%%unpack%%%%%
varName = fieldnames(opInput);
for i = 1:numel(varName)
    eval([varName{i} '=opInput.' varName{i} ';']);
end
fitFxn = @(x)leastSQ(x,CPmap,concMap,f_base,f_diff,c_base,c_diff,noT,Nout_int,...
    noEMUperm,cauchyTags,Sint,deltaT,tSampleIndex,noTsample,dataMet,fluxStoicT,...
    Nout,noExpData,dataMetConc_EXPvect,dataMetSE_EXPvect,dataMetMID_SIMvect,...
    dataMetMID_EXPvect,dataMetMID_SEvect,EMUstateStoreIS_block,EMUstateStore,...
    EMUstate,stagMap,cstag,metConcProfile_SIM,A_cell,Cmap_cell,Vmap_cell,EMUsize,...
    radData);
%%%%%%%%%%%%%%%%%%%%%

xStart = xFeas;
tStart = tic;
[xFinish,fval,exitflag]= fmincon(fitFxn,xStart,opCon.Acon,opCon.Bcon,[],[],opCon.lb,opCon.ub,[],opOptions);
tElapse = toc(tStart);