%%%clone one of optimisation instance into MC instances%%%
%%%reload data and then corrupt%%%

OFspec = OpenFLUX.OFobjSpecification(mfilename,OFspecFileName);
load(strcat([OFspec.loadFolder, OFspec.fileName]));

OF.isMonteCarlo = true;
OF.mcCloneSource = OFspec.fileName;

opSave.xFeas = [];
opSave.solnIndex = 0;
opSave.xFitSeries = struct();%datetime saved, xStart, tElapse, xFinish, fval, exitflag

fileListHPC = cell(OFspec.noCase*OFspec.noRepeatsPerCase,3);
cc = 1;
midTotLength = OF.dataMet{end,9}(end,2);
for ixx = 1:OFspec.noCase
    %%%corrupt data%%%
    while 1
        for j = 1:size(OF.dataMet,1)
            fracAb = OF.dataMet{j,2};
            fractAbSE = OF.dataMet{j,3};
            sizeAB = size(fracAb);
            %     %%%need to populate zero values before calculating fractions%%%
            fracAb(fracAb<0) = 0;
            for k = 1:sizeAB(2)
                abUse = fracAb(:,k);
                errUse = fractAbSE(:,k);
                allMat = abUse + errUse.*randn(sizeAB(1),1);
                allMat(allMat<0) = 0;
                allMat_sum = sum(allMat);
                allMat_fract = allMat/allMat_sum;
                OF.dataMet{j,4}(k) = allMat_sum;
                OF.dataMet{j,6}(:,k) = allMat_fract;
            end
        end
        
        dataMetConc_EXP = cell2mat(OF.dataMet(:,4));
        dataMetConc_EXPvect = dataMetConc_EXP(:);
        dataMetMID_EXPvect = zeros(midTotLength,1);
        for j = 1:size(OF.dataMet,1)
            dataMetMID_EXPvect(OF.dataMet{j,9}(:,2)) = OF.dataMet{j,6}(OF.dataMet{j,9}(:,1));
        end
        if ~any(isnan(dataMetMID_EXPvect))
            break
        end
    end
    %specify extra data (radiation)
    OFadditionalDataScript%placed here to corrupt it once per case
    OF.additionalData = additionalData;
    %%%%%%
    
    for j = 1:OFspec.noRepeatsPerCase
        OF.intKntPos = sort(rand(size(OF.intKntPos)));%%%guess rand knot
        simParas = OF.prepSimulation;
        x0 = rand(size(simParas.lb));
        fitFxn = OF.generateFitFxn;
        opSave.datetimeCreated = datetime('now','Format','yyyyMMdd_HHmmSSS');
        opSave.fitFxn = fitFxn;
        opSave.x0 = x0;
        opSave.xIntKnot = OF.intKntPos;
        OF.mcCaseRep = [ixx j];
        
        if OF.isODEsolver
            opSave.isODE = true;
            opSave.conFxn = simParas.conFxn;
            opSave.saveFileName = strcat(['ODEmc_' char(opSave.datetimeCreated) '.mat']);
        else
            opSave.isODE = false;
            opSave.stepBTWsample = OF.stepBTWsample;
            opSave.AconParas = simParas.AconParas;
            opSave.saveFileName = strcat(['SBRmc_' char(opSave.datetimeCreated) '.mat']);
        end
        
        fileListHPC{cc,1} = opSave.saveFileName;
        fileListHPC{cc,2} = ixx;
        fileListHPC{cc,3} = j;
        cc = cc + 1;
        save(strcat(OFspec.mcSavFolder, opSave.saveFileName), 'opSave','OF');
        disp(strcat(OFspec.mcSavFolder, opSave.saveFileName));
        pause(0.2);
    end
end
save HPCmcFileList fileListHPC