%%%specify measurement data%%%
%%%modify optimisation parameters%%%
%%%run optimisation%%%%
clearvars
clc
rng('shuffle');
%%%%load ins soln%%%%
aa = load('old/ins_HR_186.mat');


load OFobj
OF.isODEsolver = true;
OF.isOptimisation = true;
OF.genLabelledSubstrate;
OF.metDataFileName = 'metDataIns';%for insulin
% OF.metDataFileName = 'metDataBas.txt';%for basal
reLoadData = true;
% reLoadData = false;
%%%setup data input%%%
%run this for the first time to import and estimate data error
%reload same dataset to ensure optimisation is consistently parameterised  (initial conc)
%%%%this is for insulin%%%%
% %{
% %{
if reLoadData    
    OF.reEstimateError('load',[]);%%reload previous error estimates
else
    OF.importMetData;
    OF.reEstimateError('generate',10000);%%run MC to estimate error
    OF.reEstimateError('save',[]);
end

%constrain accoa_out and glycogen_out initial to near zero 0.01
hitMet = strcmp('ACCOA_out',OF.metListInt);
OF.concScale(hitMet,:) = 0.01;
hitMet = strcmp('GLYCOGEN_out',OF.metListInt);
OF.concScale(hitMet,:) = 0.01;

%%%%setup optimisation%%%%
OF.intKntPos = sort(rand(size(OF.intKntPos)));%%%guess rand knot
OF.intKntPos = aa.xKnot;
simParas = OF.prepSimulation;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%specify extra data (radiation)
%%%radiolabel data%%%
%%%setting up correction matrix
%this is for FA accumulation data, as acetyl-CoA
%c1 glucose, c2 mix, c3 unlabeled
radData{1,1} = pinv([cVectGen2([0.01 0.99],2,3)'
    cauchy(cVectGen2([1-0.0107 0.0107],1,2)',cVectGen2([0.01 0.99],1,2)')
    cVectGen2([1-0.0107 0.0107],2,3)']');
%this is for glycogen, as glucose monomer
radData{2,1} = pinv([cVectGen2([0.01 0.99],6,7) cVectGen2([0.9893 0.0107],6,7)]);

%%%specifying data and data position
%%FA accumulation pCmol/mgP accumulated for 1 hr
if strcmp('metDataIns',OF.metDataFileName)%for insulin
    radData{1,2} = [76973.028 9680.29035];%average ,se for FA
    radData{2,2} = [241.2611 10.16552]*6000;%average ,se for glycogen
else %for basal
    radData{1,2} = [16819.87	652.2888155];%average ,se for FA
    radData{2,2} = [6.413856	0.478499]*6000;%average ,se for glycogen
end
radData{1,3} = OF.findEMUindex('ACCOA_out',[1 1]);
radData{2,3} = OF.findEMUindex('GLYCOGEN_out',[1 1 1 1 1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = rand(size(simParas.lb));
% x0 = aa.xfitSeries{end,3};
% opOptions = optimoptions('fmincon', 'Display','iter','MaxFunEvals',40000);
% xGuess = x0;
% while 1
%     xFeas = fmincon(@(x)minDistX0(x,x0),xGuess,[],[],[],[],simParas.lb,simParas.ub,simParas.conFxn,opOptions);
%     if all(simParas.conFxn(xFeas)<=0)
%         break
%     else
%         xGuess = xFeas;
%     end
% end
% return

fitFxn = OF.generateFitFxn(radData);
%%%%%%%%%%%%%%%%%%%%%
xFeas = aa.xfitSeries{end,3};
disp(fitFxn(xFeas));

xStart = xFeas;
tStart = tic;
[xFinish,fval,exitflag]= fmincon(fitFxn,xStart,[],[],[],[],simParas.lb,simParas.ub,simParas.conFxn,opOptions);
tElapse = toc(tStart);

return


%%%other useful functions%%%%
%%%%unpack%%%%%
%{
varName = fieldnames(opInput);
for i = 1:numel(varName)
    eval([varName{i} '=opInput.' varName{i} ';']);
end
%}