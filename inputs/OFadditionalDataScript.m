%%%radiolabel data%%%
%%%setting up correction matrix
%this is for FA accumulation data, as acetyl-CoA
%c1 glucose, c2 mix, c3 unlabeled
additionalData{1,1} = pinv([OpenFLUX.cVectGen2([0.01 0.99],2,3)'
    cauchy(OpenFLUX.cVectGen2([1-0.0107 0.0107],1,2)',OpenFLUX.cVectGen2([0.01 0.99],1,2)')
    OpenFLUX.cVectGen2([1-0.0107 0.0107],2,3)']');
%this is for glycogen, as glucose monomer
additionalData{2,1} = pinv([OpenFLUX.cVectGen2([0.01 0.99],6,7) OpenFLUX.cVectGen2([0.9893 0.0107],6,7)]);

%%%specifying data and data position
%%FA accumulation pCmol/mgP accumulated for 1 hr
if strcmp('metDataIns',OF.metDataFileName)%for insulin
    additionalData{1,2} = [76973.028 9680.29035];%average ,se for FA
    additionalData{2,2} = [241.2611 10.16552]*6000;%average ,se for glycogen
else %for basal
    additionalData{1,2} = [16819.87	652.2888155];%average ,se for FA
    additionalData{2,2} = [6.413856	0.478499]*6000;%average ,se for glycogen
end
additionalData{1,3} = OF.findEMUindex('ACCOA_out',[1 1]);
additionalData{2,3} = OF.findEMUindex('GLYCOGEN_out',[1 1 1 1 1 1]);