%put together optimisation results from normal optimisations or MC batch
%plot distribution of optimisation results: flux, fval, knots
%check all instances in a folder for optimisation with same specifications
%and feasible solution as the reference optimisation instance

OFspec = OpenFLUX.OFobjSpecification(mfilename,OFspecFileName);
if OFspec.reloadCompiledResults
    load(OFspec.compileFileName);
else
    opRef = load(strcat([OFspec.loadFolder OFspec.opReferenceFile]));
    fileList = dir(strcat([OFspec.loadFolder '*.mat']));
    compiledResults = [];
    cc = 1;
    for i = 1:size(fileList,1)
        instLoad = load(strcat([OFspec.loadFolder fileList(i).name]));
        checkCrit = false(5,1);
        checkCrit(1) = numel(instLoad.opSave.lb) == numel(opRef.opSave.lb);
        checkCrit(2) = instLoad.opSave.solnIndex > 0;
        checkCrit(3) = instLoad.OF.isODEsolver == opRef.OF.isODEsolver;
        checkCrit(4) = instLoad.OF.isDynamic == opRef.OF.isDynamic;
        checkCrit(5) = numel(instLoad.OF.stepBTWsample)==numel(opRef.OF.stepBTWsample) && ...
            all(instLoad.OF.stepBTWsample == opRef.OF.stepBTWsample);
        checkCrit(6) = instLoad.OF.isMonteCarlo == opRef.OF.isMonteCarlo;
        checkCrit(7) = strcmp(instLoad.OF.modelCondition,opRef.OF.modelCondition);
        checkCrit(8) = strcmp(instLoad.OF.modelFileName, opRef.OF.modelFileName);
        
        if ~all(checkCrit)
            continue
        end
        compiledResults(cc).instanceName = fileList(i).name;
        compiledResults(cc).mcCase = instLoad.OF.mcCaseRep;
        compiledResults(cc).mcCloneSource = instLoad.OF.mcCloneSource;
        compiledResults(cc).fval = instLoad.opSave.xFitSeries(end).fval;
        compiledResults(cc).exitflag = instLoad.opSave.xFitSeries(end).exitflag;
        compiledResults(cc).xFinish = instLoad.opSave.xFitSeries(end).xFinish;
        if opRef.OF.isDynamic
            compiledResults(cc).xKnot = instLoad.opSave.xIntKnot;
            [simEMU,simConc,simFlux,simTime,mid_outParsed] = generatePlotData(instLoad.OF,compiledResults(cc).xFinish);
            compiledResults(cc).simFlux = simFlux;
            compiledResults(cc).simEMU = simEMU.emuFract;
            compiledResults(cc).mid_outParsed = mid_outParsed(:,[2 3 4]);
        else
            solnToVisualise = instLoad.opSave.xFitSeries(instLoad.opSave.solnIndex).xFinish;
            [f, EMUstate,flux,fMID_diff] = instLoad.opSave.fitFxn(solnToVisualise);
            compiledResults(cc).simFlux = flux;
        end
        
        
        cc = cc + 1;
        if opRef.OF.isDynamic
            if strcmp(strcat([OFspec.opReferenceFile '.mat']), fileList(i).name)
                emuList_ref = simEMU.emuList;
                simTime_ref = simTime;
                mid_outParsed_ref = mid_outParsed;
                sampleTime_ref = opRef.OF.sampleTime;
            end
        else
            emuList_ref = []; 
            simTime_ref = []; 
            mid_outParsed_ref = []; 
            sampleTime_ref = [];
        end
    end
        
    save(OFspec.compileFileName, 'compiledResults','simTime_ref','emuList_ref','mid_outParsed_ref','sampleTime_ref');
end

if OFspec.mcBestSoln
    caseKeep = false(size(compiledResults,2),1);
    mcCase = [];
    fvalMapped = [];
    for i = 1:size(compiledResults,2)
        if ~isempty(compiledResults(i).mcCase)
            mcCase(end+1,:) = compiledResults(i).mcCase;
            fvalMapped(end+1,1) = compiledResults(i).fval;
        end
    end
    uniqueCase = unique(mcCase(:,1));
    for i = 1:numel(uniqueCase)
        hitCases = find(mcCase(:,1)==uniqueCase(i));
        [minFval, index] = min(fvalMapped(hitCases));
        caseKeep(hitCases(index)) = true;        
    end
    compiledResults = compiledResults(caseKeep);
end


if ~opRef.OF.isDynamic
   caseStore = [];
   for i = 1:size(compiledResults,2)
       caseStore(i,:) = [compiledResults(i).mcCase compiledResults(i).fval];
   end
   caseUnique = unique(caseStore(:,1));
   fluxStore = [];
   minFStore = [];
   for i = 1:numel(caseUnique)
       hitRows = find(caseStore(:,1)==caseUnique(i));
       [minF index] = min(caseStore(hitRows,3));
       minFStore(1,i) = minF;
       fluxStore(:,i) = compiledResults(hitRows(index)).simFlux;
   end
    return
end

%draw fval distribution relative to knots in panels
%plotting all results instead of selecting best out of a set.
f1=figure(1);
xKnotTable = [];
fvalTable = [];
for i = 1:numel(compiledResults)
    xKnotTable(i,:) = compiledResults(i).xKnot;
    fvalTable(i) = compiledResults(i).fval;
end
for i = 1:size(xKnotTable,2)
    ax1 = subplot(1,size(xKnotTable,2),i);
    plot(ax1,xKnotTable(:,i),fvalTable,'.');
    xlim([0,1]);
    if i == 1
        ylabel('fval');
    end
    xlabel(strcat(['knot ' num2str(i)]));
end
f1.Position(1) = 120; f1.Position(2) = 507;
drawFig2(1,[],mid_outParsed_ref,compiledResults,simTime_ref,sampleTime_ref);
