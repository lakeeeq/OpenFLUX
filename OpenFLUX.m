classdef OpenFLUX < handle
    %author: Lake-Ee Quek (2018)
    %email: lake-ee.quek@sydney.edu.au
    %OpenFLUX for 13C-DMFA
    properties
        modelFileName %model text file
        ionFormFileName %met's name, size & ion formula
        metDataFileName %metabolite data file
        modelCondition %model comment
        orderS = 3 %order of b-spline
        noIntKnots = 1 %number of itnernal knots
        intKntPos %knot position(s)
        sampleTime %specify sampling time points
        stepBTWsample %step size configuration
        noSteps
        isDynamic = true %set false to output steady-state model
        isDEsolver = false %for ODE solver (not implemented)
        rxnEQfull %reactions in model file
        rxnEQ %reaction equation
        excludedMetabolites %list external species
        simulatedMDVs %list simulated EMUs
        metList %full metaboltie list
        metListInt %balanced metabolites
        EMUrxnList %EMU reactions
        natSub13Cenrich = 0.0107 %enrichment of other substrates
        natEndo13Cenrich = 0.0107 %enrichment of endogenous met
        labelledSub %specify input substrate name and positional enrichment
        ionForm %ion formula from txt file
        dataMet %metabolite tot abs and mass fract with estimated errors
        fluxBound = [10 1e5] %global min and max bound of fluxes
        fluxScale %flux min max for individual reactions
        concBound = [0.01 1e4] %global min and max bound of concentrations
        concScale %conc min max for individual balanced mets
        midMinError = 0.01 %minimum error of mass fraction
        OPoptions %optimisation options
    end
    properties (Access = private)
%         isUsingModelDynamic
%         isUsingSolverDE
        EMUmodelOutput
        EMUinputSubstrates
        EMUinsubMID
        EMUstateStoreIS_block %input sub state expanded over time
        EMUstate %instantenous EMU
        cauchyTags %EMUs that need cauchy productf
        matchExt %external species in 'metList'
        bigEMUmodel %information of emu balances
        knotSeq
        knotSeqInt
        Nout
        Nout_int
        tSampleIndex        
        simTime
        noIntMets
        noReactions %number of reactions in model
        noParaPerFlux %how many CP per flux   
        Sfull %stoic mat include external species
        Sint %balanced stoic mat      
    end
    events
    end
    methods
        function ofOBJ = OpenFLUX()
        end
        
        function buildModel(ofOBJ)
            %read and generate EMU model
            [rxnEQ, ofOBJ.excludedMetabolites, ofOBJ.simulatedMDVs, ofOBJ.rxnEQfull] = rxnExtractor(ofOBJ.modelFileName);
            ofOBJ.noReactions = size(rxnEQ,1);
            ofOBJ.fluxScale = [ofOBJ.fluxBound(1)*ones(ofOBJ.noReactions,1)	ofOBJ.fluxBound(2)*ones(ofOBJ.noReactions,1)];
            ofOBJ.rxnEQ = genRxnEQ(rxnEQ);
            [ofOBJ.Sfull, balRxn, ofOBJ.metList, ofOBJ.matchExt] = buildStoic(rxnEQ, ofOBJ.excludedMetabolites);
            ofOBJ.Sint = ofOBJ.Sfull(:,~ofOBJ.matchExt)';
            ofOBJ.metListInt = ofOBJ.metList(~ofOBJ.matchExt);
            noIntMets = numel(ofOBJ.metListInt);
            ofOBJ.noIntMets = noIntMets;
            ofOBJ.concScale = [ones(noIntMets,1)*ofOBJ.concBound(1) ones(noIntMets,1)*ofOBJ.concBound(2)];
            if find(~balRxn,1,'first')<find(balRxn,1,'last')
                disp('put all flux balance rxn first and S-type rxn last');
                return
            end
            [EMUrxnList, traceRxn] = buildEMUrxn(rxnEQ);
            EMUrxnPicked = pickEMUrxn(EMUrxnList, ofOBJ.simulatedMDVs);
            ofOBJ.EMUrxnList = EMUrxnList(EMUrxnPicked,:);
            ofOBJ.EMUinputSubstrates = getEMUinputSubstrates(EMUrxnList(EMUrxnPicked,:));
            if ofOBJ.isDynamic
                ofOBJ.EMUmodelOutput = buildEMUdnsMat(EMUrxnList(EMUrxnPicked,:),ofOBJ.EMUinputSubstrates,...
                    ofOBJ.Sfull, ofOBJ.metList,ofOBJ.matchExt);
                %print SS model
            else
                [EMUbalanceBlock, EMUcalculated, EMUsimulated_out, EMUrxnList_collapsed] =...
                    buildEMUssMat(EMUrxnList(EMUrxnPicked,:), simulatedMDVs, EMUinputSubstrates);
                printSSmfiles(EMUbalanceBlock, EMUsimulated_out, EMUinputSubstrates);
            end
        end
        
        function genLabelledSubstrate(ofOBJ,substrateLabelled)
            %build input substrate MIDs
            ofOBJ.labelledSub = substrateLabelled;
            ofOBJ.EMUinsubMID = inputSubBuilder(ofOBJ);
        end
        
        function importMetData(ofOBJ)
            %read met data file
            ofOBJ.dataMet = readMetDatFile(ofOBJ.metDataFileName);
        end
        
        function reEstimateError(ofOBJ,actionToDo,noItt,matFileName)
            %actionToDo: 'load', 'save', 'generate'
            switch actionToDo
                case 'load'
                    load(matFileName,'-mat');
                    ofOBJ.dataMet = dataMet;
                case 'save'
                    dataMet = ofOBJ.dataMet;
                    save(matFileName,'dataMet');
                case 'generate'
                    ofOBJ.dataMet = calcErrorByMC(ofOBJ.dataMet,noItt);
            end
            genConcScale(ofOBJ);%%%specify conc range
            
        end
        
        function [opCon,opInput] = prepSimulation(ofOBJ)
            [ofOBJ.noSteps, ofOBJ.simTime, ofOBJ.tSampleIndex, ofOBJ.noParaPerFlux,...
                ofOBJ.knotSeq, ofOBJ.knotSeqInt, ofOBJ.Nout, ofOBJ.Nout_int] = genBSplineMat(ofOBJ);
            
            [ofOBJ.bigEMUmodel, ofOBJ.cauchyTags, ofOBJ.EMUstateStoreIS_block, ofOBJ.EMUstate]...
                = genEMUmodelStart(ofOBJ);
            ofOBJ.ionForm = readIonFormFile(ofOBJ.ionFormFileName);
            ofOBJ.dataMet = expandDataMet(ofOBJ);
            [opCon,opInput] = genOptimisationProblem(ofOBJ);
        end
                
        function corruptData()
        end
        
        function varOut = getPrivProp(ofOBJ,varName)
            %shortcut to access private properties from outside
            varOut = eval(['ofOBJ.',varName]);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rxnEQ, excludedMetabolites, simulatedMDVs,rxnEQfull] = rxnExtractor(fileName)
%read model from file
rxnEQ = '';
rxnEQfull = '';
excludedMetabolites = {};
simulatedMDVs = {};
fid = fopen(fileName,'r');
modelFile = textscan(fid,'%s','delimiter','\n');
modelFile = modelFile{1};
fclose(fid);

rc = 1;
while rc <= size(modelFile,1)
    if isempty(modelFile{rc})
        break
    end
    rc = rc+1;
end
rxnText = modelFile(2:rc-1);
rxnEQ = cell(size(rxnText,1),4);
for i = 1:size(rxnText,1)
    [coeff, met, ctrans, clabel, rxnBasis] = line2species(rxnText{i});
    rxnEQ{i,1} = coeff;
    rxnEQ{i,2} = met;
    rxnEQ{i,3} = ctrans;
    rxnEQ{i,4} = clabel;
    rxnEQ{i,5} = rxnBasis;
end


while rc <= size(modelFile,1)
    if ~isempty(modelFile{rc}) && all(modelFile{rc}(1:2)=='##')
        break
    else
        rc = rc + 1;
    end
end
listType = modelFile{rc};
rc = rc + 1;
metList = {};
while rc <= size(modelFile,1)
    metLine = modelFile{rc};
    if isempty(regexp(metLine,'#\t','ONCE'))
        break
    end
    tabPos = regexp(metLine,'\t');
    if numel(tabPos)== 1
        metList{end+1,1} = metLine(tabPos+1:end);
    else
        metList{end+1,1} = metLine(tabPos(1)+1:tabPos(2)-1);
    end
    rc = rc + 1;
end
metList = strtrim(metList);
if ~isempty(regexpi(listType,'excludedMetabolites'))
    excludedMetabolites = metList;
elseif ~isempty(regexpi(listType,'simulatedMDVs'))
    simulatedMDVs = metList;
end


while rc <= size(modelFile,1)
    if ~isempty(modelFile{rc}) && all(modelFile{rc}(1:2)=='##')
        break
    else
        rc = rc + 1;
    end
end
listType = modelFile{rc};
rc = rc + 1;
metList = {};
while rc <= size(modelFile,1)
    metLine = modelFile{rc};
    if isempty(regexp(metLine,'#\t','ONCE'))
        break
    end
    tabPos = regexp(metLine,'\t');
    if numel(tabPos)== 1
        metList{end+1,1} = metLine(tabPos+1:end);
    else
        metList{end+1,1} = metLine(tabPos(1)+1:tabPos(2)-1);
    end
    rc = rc + 1;
end
metList = strtrim(metList);
if ~isempty(regexpi(listType,'excludedMetabolites'))
    excludedMetabolites = metList;
elseif ~isempty(regexpi(listType,'simulatedMDVs'))
    simulatedMDVs = metList;
end


%test print reaction network
for i = 1:size(rxnEQ,1)
    hitProd = rxnEQ{i,1}>0;
    hitReact = rxnEQ{i,1}<0;
    fprintf('R%g\t',i);
    LHS = '';RHS = '';
    for j = find(hitProd)
        if isempty(RHS)
            RHS = strcat([num2str(rxnEQ{i,1}(j)), ' ', rxnEQ{i,2}{j}, '(',...
                rxnEQ{i,3}{j},')']);
        else
            RHS = strcat([RHS,' + ', num2str(rxnEQ{i,1}(j)), ' ', rxnEQ{i,2}{j},...
                '(',rxnEQ{i,3}{j},')']);
        end
    end
    for j = find(hitReact)
        if isempty(LHS)
            LHS = strcat([num2str(-rxnEQ{i,1}(j)), ' ', rxnEQ{i,2}{j}, '(',...
                rxnEQ{i,3}{j},')']);
        else
            LHS = strcat([LHS,' + ', num2str(-rxnEQ{i,1}(j)), ' ', rxnEQ{i,2}{j},...
                '(',rxnEQ{i,3}{j},')']);
        end
    end
    fprintf('%s = %s\n',LHS,RHS);
    rxnEQfull{end+1,1} = strcat([LHS ' = ' RHS]);
end
fprintf('\n');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coeff, met, ctrans, rxnType, rxnBasis] = line2species(rxnLine)
%extract OpenFLUX model reactions, metabolite and atom transition line by
%line
coeff = [];
met = {};
ctrans = {};
rxnType = '';
rxnBasis = '';

tabPos = regexp(rxnLine,'\t');
rxnEQline_r = strtrim(rxnLine(tabPos(1)+1:tabPos(2)-1));
rxnEQline_t = strtrim(rxnLine(tabPos(2)+1:tabPos(3)-1));

%in case line doesn't end with reaction type, take care of extra tabs
if numel(tabPos) <= 4
    rxnType = rxnLine(tabPos(4)+1:end);
else
    rxnType = rxnLine(tabPos(4)+1:tabPos(5)-1);
end
if numel(tabPos) <= 5
    rxnBasis = rxnLine(tabPos(5)+1:end);
else
    rxnBasis = rxnLine(tabPos(5)+1:tabPos(6)-1);
end
hitEqual = regexp(rxnEQline_r,'=','ONCE');

%do this for LHS of the equation, reaction
%negative stoic
rxnEQline = strtrim(rxnEQline_r(1:hitEqual-1));
while ~isempty(rxnEQline)
    hitPos = regexp(rxnEQline,'+','ONCE');
    %this lot do not have + sign anymore
    if isempty(hitPos)
        rxnEQpart = strtrim(rxnEQline);
        hitSpace = regexp(rxnEQpart,' ','ONCE');
        if ~isempty(hitSpace)
            stoic = strtrim(rxnEQpart(1:hitSpace-1));
            coeff(end+1) = -str2num(stoic);
            met{end+1} = strtrim(rxnEQpart(hitSpace+1:end));
        else
            coeff(end+1) = -1;
            met{end+1} = strtrim(rxnEQpart);
        end
        break
    end
    
    %metabolite and stoic must be 1 space apart
    rxnEQpart = strtrim(rxnEQline(1:hitPos-1));
    if any(rxnEQpart==' ')
        hitSpace = regexp(rxnEQpart,' ','ONCE');
        stoic = strtrim(rxnEQpart(1:hitSpace-1));
        coeff(end+1) = -str2num(stoic);
        met{end+1} = strtrim(rxnEQpart(hitSpace+1:end));
    else
        coeff(end+1) = -1;
        met{end+1} = strtrim(rxnEQpart);
    end
    rxnEQline(1:hitPos) = [];
end



%do this for RHS of the equation, reaction
%positive stoic
rxnEQline = strtrim(rxnEQline_r(hitEqual+1:end));
while ~isempty(rxnEQline)
    hitPos = regexp(rxnEQline,'+','ONCE');
    %this lot do not have + sign anymore
    if isempty(hitPos)
        rxnEQpart = strtrim(rxnEQline);
        hitSpace = regexp(rxnEQpart,' ','ONCE');
        if ~isempty(hitSpace)
            stoic = strtrim(rxnEQpart(1:hitSpace-1));
            coeff(end+1) = str2num(stoic);
            met{end+1} = strtrim(rxnEQpart(hitSpace+1:end));
        else
            coeff(end+1) = 1;
            met{end+1} = strtrim(rxnEQpart);
        end
        break
    end
    
    %metabolite and stoic must be 1 space apart
    rxnEQpart = strtrim(rxnEQline(1:hitPos-1));
    if any(rxnEQpart==' ')
        hitSpace = regexp(rxnEQpart,' ','ONCE');
        stoic = strtrim(rxnEQpart(1:hitSpace-1));
        coeff(end+1) = str2num(stoic);
        met{end+1} = strtrim(rxnEQpart(hitSpace+1:end));
    else
        coeff(end+1) = 1;
        met{end+1} = strtrim(rxnEQpart);
    end
    rxnEQline(1:hitPos) = [];
end

%%%%%%%%%
%now do for atom transition equation, if available
if isempty(rxnEQline_t)
    ctrans = met;
    for i = 1:numel(ctrans)
        ctrans{i} = '';
    end
    return
end


hitEqual = regexp(rxnEQline_t,'=','ONCE');

%do this for LHS of the equation, reaction
%negative stoic
rxnEQline = strtrim(rxnEQline_t(1:hitEqual-1));
while ~isempty(rxnEQline)
    hitPos = regexp(rxnEQline,'+','ONCE');
    %this lot do not have + sign anymore
    if isempty(hitPos)
        rxnEQpart = strtrim(rxnEQline);
        hitSpace = regexp(rxnEQpart,' ','ONCE');
        if ~isempty(hitSpace)
            ctrans{end+1} = strtrim(rxnEQpart(hitSpace+1:end));
        else
            ctrans{end+1} = strtrim(rxnEQpart);
        end
        break
    end
    
    %metabolite and stoic must be 1 space apart
    rxnEQpart = strtrim(rxnEQline(1:hitPos-1));
    if any(rxnEQpart==' ')
        hitSpace = regexp(rxnEQpart,' ','ONCE');
        ctrans{end+1} = strtrim(rxnEQpart(hitSpace+1:end));
    else
        ctrans{end+1} = strtrim(rxnEQpart);
    end
    rxnEQline(1:hitPos) = [];
end



%do this for RHS of the equation, reaction
%positive stoic
rxnEQline = strtrim(rxnEQline_t(hitEqual+1:end));
while ~isempty(rxnEQline)
    hitPos = regexp(rxnEQline,'+','ONCE');
    %this lot do not have + sign anymore
    if isempty(hitPos)
        rxnEQpart = strtrim(rxnEQline);
        hitSpace = regexp(rxnEQpart,' ','ONCE');
        if ~isempty(hitSpace)
            ctrans{end+1} = strtrim(rxnEQpart(hitSpace+1:end));
        else
            ctrans{end+1} = strtrim(rxnEQpart);
        end
        break
    end
    
    %metabolite and stoic must be 1 space apart
    rxnEQpart = strtrim(rxnEQline(1:hitPos-1));
    if any(rxnEQpart==' ')
        hitSpace = regexp(rxnEQpart,' ','ONCE');
        ctrans{end+1} = strtrim(rxnEQpart(hitSpace+1:end));
    else
        ctrans{end+1} = strtrim(rxnEQpart);
    end
    rxnEQline(1:hitPos) = [];
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [S, balRxn, metList, matchExt] = buildStoic(rxnEQ, excludedMetabolites)

%pick out balanceable reactions
balRxn = true(size(rxnEQ,1),1);
for i = 1:size(rxnEQ,1)
    if ~isempty(regexpi(rxnEQ{i,4},'S'))
        balRxn(i) = false;
    end
end
rxnEQbal = rxnEQ(balRxn,:);

%construct metabolite list
metList = {};
for i = 1:size(rxnEQbal,1)
    for j = 1:numel(rxnEQbal{i,2})
        metList{end+1,1} = rxnEQbal{i,2}{j};
    end
    if ~isempty(regexpi(rxnEQbal{i,4},'s'))
        balRxn(i) = false;
    end
end
metList = unique(metList);

%tag external metabolites
matchExt = false(size(metList));
extAbsent = false(size(excludedMetabolites));
for i = 1:numel(excludedMetabolites)
    hit = strcmp(excludedMetabolites{i},metList);
    if any(hit)
        matchExt(hit) = true;
    else
        extAbsent(i) = true;
    end
end
if any(extAbsent)
    fprintf('excludedMetabolites missing from network:\n');
    for i = find(extAbsent)'
        fprintf('%s\n',excludedMetabolites{i});
    end
    fprintf('\n');
end

%construct stoic matrix
%reaction rows, metabolite columns
S = zeros(size(rxnEQbal,1),size(metList,1));
for i = 1:size(rxnEQbal,1)
    for j = 1:numel(rxnEQbal{i,2})
        hitMet = strcmp(rxnEQbal{i,2}{j},metList);
        S(i,hitMet) = S(i,hitMet) + rxnEQbal{i,1}(j);
    end
end

hitOneMet = sum(S~=0,1)==1;
fprintf('singleton metabolites in S, should match excluded metabolites\n');
for i = find(hitOneMet)
    fprintf('%s\n',metList{i});
end
fprintf('\n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EMUbalanceBlock, EMUcalculated, EMUsimulated_out, EMUrxnList_collapsed]= buildEMUssMat(EMUrxnList, simulatedMDVs, EMUinputSubstrates)
%build steady-state EMU model/matrix
% clearvars
% clc
% load saveDat
% EMUrxnList = EMUrxnList(EMUrxnPicked,:);

%build simulated EMU list
EMUsimulated = cell(0,3);
for i = 1:numel(simulatedMDVs)
    simMDV = strtrim(simulatedMDVs{i});
    hashPos = regexp(simMDV,'#','ONCE');
    EMUsimulated{i,1} = simMDV(1:hashPos-1);
    EMUsimulated{i,2} = +(simMDV(hashPos+1:end)=='1');
    EMUsimulated{i,3} = numel(EMUsimulated{i,2});
end
EMUsimulated_out = EMUsimulated;

EMUrxnList_emuSize = zeros(size(EMUrxnList,1),1);
for i = 1:size(EMUrxnList,1)
    EMUrxnList_emuSize(i) = sum(EMUrxnList{i,5});
end


%%%%%%%%%%%%%
%%%collapse EMUrxnList for 1-to-1 mappings
%do this iteratively
cc = 1;
while cc <= size(EMUrxnList,1)
    EMUPname = EMUrxnList{cc,4};
    EMUPtag = EMUrxnList{cc,5};
    hitName = strcmp(EMUPname,EMUsimulated(:,1));
    
    if any(matchEMU(EMUPname,EMUPtag,EMUsimulated))
        cc = cc + 1;
        continue
    end
    
    hitRow = matchEMU(EMUPname,EMUPtag,EMUrxnList(:,[4,5]));
    EMUr = {};
    EMUr_rxn = [];
    %collect all possible reactants
    for j = find(hitRow)'
        EMUr_rxn(end+1,1) = j;
        for k = 1:numel(EMUrxnList{j,6})
            EMUr{end+1,1} = EMUrxnList{j,6}{k};
            EMUr{end,2} = EMUrxnList{j,7}{k};
            EMUr{end,3} = numel(EMUr{end,2});
        end
    end
    
    %check if only one reactant hit
    EMUr = removeEMUduplicate(EMUr);
    if size(EMUr,1)~=1
        cc = cc + 1;
        continue
    end
    
    %now find all matching reactants of product EMU, and replace product with EMUr
    %eg, if GLC match uniquely to GLC_ext, then all "GLC" reactant EMU
    %replaced with "GLC_ext"
    for j = 1:size(EMUrxnList,1)
        for k = 1:numel(EMUrxnList{j,6})
            if strcmp(EMUrxnList{j,6}{k}, EMUPname) && all(EMUrxnList{j,7}{k}==EMUPtag)
                EMUrxnList{j,6}{k} = EMUr{1};
                EMUrxnList{j,7}{k} = EMUr{2};
            end
        end
    end
    %then delete reaction
    EMUrxnList(EMUr_rxn,:) = [];
    EMUrxnList_emuSize(EMUr_rxn) = [];
end

%remove identical reactant product, no info
cc = 1;
while cc <= size(EMUrxnList,1)
    if numel(EMUrxnList{cc,6}) == 1
        if strcmp(EMUrxnList{cc,4},EMUrxnList{cc,6})
            if all(EMUrxnList{cc,5}==EMUrxnList{cc,7}{1})
                EMUrxnList(cc,:) = [];
                EMUrxnList_emuSize(cc) = [];
                continue
            end
        end
    end
    cc = cc + 1;
end

EMUrxnList_collapsed = EMUrxnList;

%%%%%%%%%%%%%
%sort by size, largest tracked first
%%put first largest EMUsimulated into EMUtoTrack when empty
%%track (construct matrix) reactant matched to EMU
%drop if match EMUcalculated
%put into EMUcalc if unknown of same size
%put into EMUknown if known substrate or smaller size
%put these reactants into EMU to track
%and remove these reactants from EMUsimulated if match

%size, A, X, empty A, B, Y, empty B (AX = BY), x is calculated
EMUbalanceBlock = cell(0,7);
EMUcalculated = cell(0,3);
while ~isempty(EMUsimulated)
    [maxEMU, hitEMU] = max(cell2mat(EMUsimulated(:,3)));
    EMUtoTrack = EMUsimulated(hitEMU(1),:);
    EMUsimulated(hitEMU(1),:) = [];
    while ~isempty(EMUtoTrack)
        maxEMUstart = max(cell2mat(EMUtoTrack(:,3)));
        
        X_block = {}; %product/unknown block
        Y_block = {}; %reactant/known block
        EMUcalculatedNew = cell(0,3);
        pRow = 0; rCol = 0;
        %build known/unknown variables X, Y
        while 1 %keep doing EMUs of the same size, take first
            EMUtoTrack = removeEMUduplicate(EMUtoTrack);
            hitEMU_tt = find(cell2mat(EMUtoTrack(:,3))==maxEMUstart,1,'first');
            if isempty(hitEMU_tt)
                break
            end
            EMUname = EMUtoTrack{hitEMU_tt,1};
            EMUtag = EMUtoTrack{hitEMU_tt,2};
            pRow = pRow + 1;
            X_block{pRow,1} = EMUname;
            X_block{pRow,2} = EMUtag;
            
            if ~any(matchEMU(EMUname,EMUtag,EMUcalculated))
                EMUcalculatedNew(end+1,:) = EMUtoTrack(hitEMU_tt,:);
            end
            
            %check if repeat of emu simulated, then delete
            hitRepeat = matchEMU(EMUname,EMUtag,EMUsimulated);
            if any(hitRepeat)
                EMUsimulated(hitRepeat,:) = [];
            end
            
            EMUtoTrack(hitEMU_tt,:) = [];
            hitEMUs_p = matchEMU(EMUname,EMUtag,EMUrxnList(:,[4 5]));
            for i = find(hitEMUs_p)'%matched reaction with such product
                %now check reactants
                if numel(EMUrxnList{i,6}) == 1
                    %match known
                    EMUrName = EMUrxnList{i,6}{1};
                    EMUrTag = EMUrxnList{i,7}{1};
                    if any(matchEMU(EMUrName,EMUrTag,EMUcalculated))
                        rCol = rCol + 1;
                        Y_block{rCol,1} = EMUrName;
                        Y_block{rCol,2} = EMUrTag;
                        continue
                    end
                    %match input substrate
                    if any(matchEMU(EMUrName,EMUrTag,EMUinputSubstrates))
                        rCol = rCol + 1;
                        Y_block{rCol,1} = EMUrName;
                        Y_block{rCol,2} = EMUrTag;
                        continue
                    end
                    %else if unknown
                    pRow = pRow + 1;
                    X_block{pRow,1} = EMUrName;
                    X_block{pRow,2} = EMUrTag;
                    
                    if ~any(matchEMU(EMUrName,EMUrTag,EMUcalculatedNew))
                        EMUtoTrack{end+1,1} = EMUrName;
                        EMUtoTrack{end,2} = EMUrTag;
                        EMUtoTrack{end,3} = sum(EMUrTag);
                    end
                else
                    %match smaller, must be a pair
                    EMUrName = EMUrxnList{i,6}{1};
                    EMUrTag = EMUrxnList{i,7}{1};
                    rCol = rCol + 1;
                    Y_block{rCol,1} = EMUrName;
                    Y_block{rCol,2} = EMUrTag;
                    hit1 = matchEMU(EMUrName,EMUrTag,EMUcalculated);
                    hit2 = matchEMU(EMUrName,EMUrTag,EMUinputSubstrates);
                    hit3 = matchEMU(EMUrName,EMUrTag,EMUcalculatedNew);
                    if ~any(hit1) && ~any(hit2) && ~any(hit3)
                        EMUtoTrack{end+1,1} = EMUrName;
                        EMUtoTrack{end,2} = EMUrTag;
                        EMUtoTrack{end,3} = sum(EMUrTag);
                    end
                    
                    EMUrName = EMUrxnList{i,6}{2};
                    EMUrTag = EMUrxnList{i,7}{2};
                    Y_block{rCol,3} = EMUrName;
                    Y_block{rCol,4} = EMUrTag;
                    hit1 = matchEMU(EMUrName,EMUrTag,EMUcalculated(:,[1 2]));
                    hit2 = matchEMU(EMUrName,EMUrTag,EMUinputSubstrates(:,[1 2]));
                    hit3 = matchEMU(EMUrName,EMUrTag,EMUcalculatedNew);
                    if ~any(hit1) && ~any(hit2) && ~any(hit3)
                        EMUtoTrack{end+1,1} = EMUrName;
                        EMUtoTrack{end,2} = EMUrTag;
                        EMUtoTrack{end,3} = sum(EMUrTag);
                    end
                end
            end
        end
        
        %build balance matrix A, B (AX = BY)
        %remove duplicates
        X_block = removeEMUduplicate(X_block);
        if size(Y_block,2)>2
            cauchyRows = false(size(Y_block,1),1);
            for i = 1:size(Y_block,1)
                if ~isempty(Y_block{i,3})
                    cauchyRows(i) = true;
                end
            end
            Y_blockSingles = removeEMUduplicate(Y_block(~cauchyRows,:));
            Y_blockCauchy = removeEMUduplicatePair(Y_block(cauchyRows,:));
            Y_block = [Y_blockSingles;Y_blockCauchy];
        else
            Y_block = removeEMUduplicate(Y_block);
        end
        
        
        %A_block and B_block in row, column, stoic, reactionIndex
        A_block = zeros(0,4); B_block = zeros(0,4);
        A_block_size = [size(X_block,1) size(X_block,1)];
        B_block_size = [size(X_block,1) size(Y_block,1)];
        
        for i = 1:size(X_block,1)
            hitEMU = matchEMU(X_block{i,1},X_block{i,2},EMUrxnList(:,[4,5]));
            for j = find(hitEMU)'
                %match to Y_block, cauchy
                if numel(EMUrxnList{j,6})==2
                    name1 = EMUrxnList{j,6}{1}; name2 = EMUrxnList{j,6}{2};
                    tag1 = EMUrxnList{j,7}{1}; tag2 = EMUrxnList{j,7}{2};
                    for k = 1:size(Y_block,1)
                        if isempty(Y_block{k,3})
                            continue
                        end
                        testname1 = Y_block{k,1}; testname2 = Y_block{k,3};
                        testtag1 = Y_block{k,2}; testtag2 = Y_block{k,4};
                        hitYRow = false;
                        if strcmp(name1,testname1) && all(tag1==testtag1) &&...
                                strcmp(name2,testname2) && all(tag2==testtag2)
                            hitYRow = true;
                        elseif strcmp(name1,testname2) && all(tag1==testtag2) &&...
                                strcmp(name2,testname1) && all(tag2==testtag1)
                            hitYRow = true;
                        end
                        if hitYRow
                            B_block(end+1,:) = [i k -EMUrxnList{j,3} EMUrxnList{j,2}];
                        end
                    end
                    continue
                end
                %%%%%%
                
                EMUrName = EMUrxnList{j,6}{1}; EMUrTag = EMUrxnList{j,7}{1};
                
                %match to Y_block, input substrate
                if any(matchEMU(EMUrName,EMUrTag,EMUinputSubstrates(:,[1 2])))
                    for k = 1:size(Y_block,1)
                        if strcmp(EMUrName,Y_block{k,1}) && all(EMUrTag==Y_block{k,2})
                            B_block(end+1,:) = [i k -EMUrxnList{j,3} EMUrxnList{j,2}];
                        end
                    end
                    continue
                end
                %%%%%
                %match to Y_block, known
                if any(matchEMU(EMUrName,EMUrTag,EMUcalculated(:,[1 2])))
                    for k = 1:size(Y_block,1)
                        if strcmp(EMUrName,Y_block{k,1}) && all(EMUrTag==Y_block{k,2})
                            B_block(end+1,:) = [i k -EMUrxnList{j,3} EMUrxnList{j,2}];
                        end
                    end
                    continue
                end
                %%%%%
                %match to X_block, unknown
                
                for k = 1:size(X_block,1)
                    if k == i
                        continue
                    end
                    if strcmp(EMUrName,X_block{k,1}) && all(EMUrTag==X_block{k,2})
                        A_block(end+1,:) = [i k EMUrxnList{j,3} EMUrxnList{j,2}];
                    end
                end
            end
        end
        
        %add new EMUcalculated to EMUcalculated
        EMUcalculatedNew = removeEMUduplicate(EMUcalculatedNew);
        EMUcalculated = [EMUcalculated; EMUcalculatedNew];
        
        
        
        %combine stoic in X_block and Y_block, and sums X_block diagonal
        A_cell = cell(0,3);
        B_cell = cell(0,3);
        %construct diagonal sum
        for i = 1:size(X_block,1)
            hitRows_B = B_block(:,1) == i;
            hitRows_A = A_block(:,1) == i;
            v_unique = unique([B_block(hitRows_B,4);A_block(hitRows_A,4)]);
            stoicCoeff = zeros(numel(v_unique),1);
            for j = find(hitRows_B)'
                hitV = v_unique==B_block(j,4);
                stoicCoeff(hitV) = stoicCoeff(hitV) + B_block(j,3);
            end
            for j = find(hitRows_A)'
                hitV = v_unique==A_block(j,4);
                stoicCoeff(hitV) = stoicCoeff(hitV) - A_block(j,3);
            end
            A_cell{end+1,1} = [i i];
            A_cell{end,2} = stoicCoeff;
            A_cell{end,3} = v_unique;
        end
        
        %construct A_cell non diagonal
        if ~isempty(A_block)
            [cellUnique, IA, IC] = unique(A_block(:,[1,2]),'rows');
            for i = 1:size(cellUnique,1)
                hitCU = IC == i;
                v_unique = unique(A_block(hitCU,4));
                stoicCoeff = zeros(numel(v_unique),1);
                for j = find(hitCU)'
                    hitV = v_unique==A_block(j,4);
                    stoicCoeff(hitV) = stoicCoeff(hitV) + A_block(j,3);
                end
                A_cell{end+1,1} = cellUnique(i,:);
                A_cell{end,2} = stoicCoeff;
                A_cell{end,3} = v_unique;
            end
        end
        
        %construct B_cell
        if ~isempty(B_block)
            [cellUnique, IA, IC] = unique(B_block(:,[1,2]),'rows');
            
            for i = 1:size(cellUnique,1)
                hitCU = IC == i;
                v_unique = unique(B_block(hitCU,4));
                stoicCoeff = zeros(numel(v_unique),1);
                for j = find(hitCU)'
                    hitV = v_unique==B_block(j,4);
                    stoicCoeff(hitV) = stoicCoeff(hitV) + B_block(j,3);
                end
                B_cell{end+1,1} = cellUnique(i,:);
                B_cell{end,2} = stoicCoeff;
                B_cell{end,3} = v_unique;
            end
        end
        
        %size, A, X, empty A, B, Y, empty B (AX = BY), x is calculated
        EMUbalanceBlock{end+1,1} = maxEMUstart;
        EMUbalanceBlock{end,2} = A_cell;
        EMUbalanceBlock{end,3} = X_block;
        EMUbalanceBlock{end,4} = A_block_size;
        EMUbalanceBlock{end,5} = B_cell;
        EMUbalanceBlock{end,6} = Y_block;
        EMUbalanceBlock{end,7} = B_block_size;
    end
end

EMUsizes = cell2mat(EMUcalculated(:,3));
[~,index] = sort(EMUsizes);
EMUcalculated =  EMUcalculated(index,:);
EMUnames = EMUcalculated(:,1);
[~,index] = sort(EMUnames);
EMUcalculated = EMUcalculated(index,:);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printSSmfiles(EMUbalanceBlock, EMUsimulated_out, EMUinputSubstrates)
% clearvars
% clc
% load saveDat

fid = fopen('EMUModel_substrate.m','w');
fprintf(fid,'%%specify EMU input data here\n');
for i = 1:size(EMUinputSubstrates,1)
    EMUname = EMUinputSubstrates{i,1};
    EMUtag = char(EMUinputSubstrates{i,2}+'0');
    fprintf(fid,'%s_%s = [];\n',EMUname,EMUtag);
end
fclose(fid);


fid = fopen('EMUModel_simulated.m','w');
fprintf(fid,'%%simulated measurement vector\n');
fprintf(fid,'x_calc = [\n');
for i = 1:size(EMUsimulated_out,1)
    EMUname = EMUsimulated_out{i,1};
    EMUtag = char(EMUsimulated_out{i,2}+'0');
    fprintf(fid,'%s_%s''\n',EMUname,EMUtag);
end
fprintf(fid,'];\n');
fclose(fid);


fid1 = fopen('EMUModel_calculated.m','w');
fid2 = fopen('EMUModel_loader.m','w');
EMUsizes = cell2mat(EMUbalanceBlock(:,1));
EMUsize_unique = unique(EMUsizes);
for i = 1:numel(EMUsize_unique)
    hitEMUsize = EMUsizes==EMUsize_unique(i);
    rc = 1;
    for j = find(hitEMUsize)'
        Astring = strcat(['A' num2str(EMUsize_unique(i)) '_' num2str(rc)]);
        Bstring = strcat(['B' num2str(EMUsize_unique(i)) '_' num2str(rc)]);
        Ystring = strcat(['Y' num2str(EMUsize_unique(i)) '_' num2str(rc)]);
        Xstring = strcat(['X' num2str(EMUsize_unique(i)) '_' num2str(rc)]);
        AinvString = strcat(['A' num2str(EMUsize_unique(i)) '_' num2str(rc) '_inv']);
        
        EMUknown = EMUbalanceBlock{j,6};
        EMUcalc = EMUbalanceBlock{j,3};
        BblockSize = EMUbalanceBlock{j,7};
        Bblock = EMUbalanceBlock{j,5};
        AblockSize = EMUbalanceBlock{j,4};
        Ablock = EMUbalanceBlock{j,2};
        
        %print B matrix
        fprintf(fid2,'%s = [',Bstring);
        matBlock = buildBlockContent(BblockSize, Bblock);
        for m = 1:BblockSize(1)
            for n = 1:BblockSize(2)
                if n == 1
                    if isempty(matBlock{m,n})
                        fprintf(fid2,'0');
                    else
                        fprintf(fid2,'%s',matBlock{m,n});
                    end
                else
                    if isempty(matBlock{m,n})
                        fprintf(fid2,',0');
                    else
                        fprintf(fid2,',%s',matBlock{m,n});
                    end
                end
            end
            if m<BblockSize(1)
                fprintf(fid2,'\n');
            end
        end
        fprintf(fid2,'];\n');
        
        %print A matrix
        fprintf(fid2,'%s = [',Astring);
        matBlock = buildBlockContent(AblockSize, Ablock);
        for m = 1:AblockSize(1)
            for n = 1:AblockSize(2)
                if n == 1
                    if isempty(matBlock{m,n})
                        fprintf(fid2,'0');
                    else
                        fprintf(fid2,'%s',matBlock{m,n});
                    end
                else
                    if isempty(matBlock{m,n})
                        fprintf(fid2,',0');
                    else
                        fprintf(fid2,',%s',matBlock{m,n});
                    end
                end
            end
            if m<AblockSize(1)
                fprintf(fid2,'\n');
            end
        end
        fprintf(fid2,'];\n');
        
        %print inversion equation
        fprintf(fid2,'%s = %s\\%s;\n',AinvString, Astring, Bstring);
        fprintf(fid2,'\n\n');
        rc = rc + 1;
        
        
        fprintf(fid1,'%s = [',Ystring);
        if size(EMUknown,2)<=2
            for k = 1:size(EMUknown,1)
                EMUname = EMUknown{k,1};
                EMUtag = char(EMUknown{k,2}+'0');
                fprintf(fid1,'%s_%s',EMUname,EMUtag);
                if k<size(EMUknown,1)
                    fprintf(fid1,'\n');
                end
            end
        else
            for k = 1:size(EMUknown,1)
                if isempty(EMUknown{k,3})
                    EMUname = EMUknown{k,1};
                    EMUtag = char(EMUknown{k,2}+'0');
                    fprintf(fid1,'%s_%s',EMUname,EMUtag);
                else
                    EMUname1 = EMUknown{k,1};
                    EMUtag1 = char(EMUknown{k,2}+'0');
                    EMUname2 = EMUknown{k,3};
                    EMUtag2 = char(EMUknown{k,4}+'0');
                    fprintf(fid1,'cauchy(%s_%s,%s_%s)',EMUname1,EMUtag1,EMUname2,EMUtag2);
                end
                if k<size(EMUknown,1)
                    fprintf(fid1,'\n');
                end
            end
            
        end
        fprintf(fid1,'];\n');
        
        fprintf(fid1,'%s = %s*%s;\n',Xstring,AinvString,Ystring);
        
        for k = 1:size(EMUcalc,1)
            EMUname = EMUcalc{k,1};
            EMUtag = char(EMUcalc{k,2}+'0');
            fprintf(fid1,'%s_%s = %s(%s,:);\n',EMUname,EMUtag,Xstring,num2str(k));
        end
        fprintf(fid1,'\n');
    end
end
fclose(fid1);
fclose(fid2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EMUrxnList, traceRxn] = buildEMUrxn(rxnEQ)
% clearvars
% clc
% load saveDat
%rxnEQ full, rxnID, prodStoic, prodN, prodAT, reactN(rows), reactAT(rows)
EMUrxnList = cell(0,7);
emuCC = 1;
traceRxn = true(size(rxnEQ,1),1);
for i = 1:size(rxnEQ,1)
    if ~isempty(regexp(rxnEQ{i,4},'B','ONCE'))
        traceRxn(i) = false;
    end
end
rxnEQTrace = rxnEQ(traceRxn,:);

%emu rxn, map product to rxn
%construct full carbon atom interaction map
%rxn#, product name, reactant name, product AT, reactant AT, product mappings, reactant mappings
metIntMap = cell(0,7);
rc = 1;
for i = 1:size(rxnEQTrace,1)
    hitProd = find(rxnEQTrace{i,1}>0);
    hitReact = find(rxnEQTrace{i,1}<0);
    for j = hitProd
        prodN = rxnEQTrace{i,2}{j};
        prodAT = rxnEQTrace{i,3}{j};
        if ~isempty(regexp(prodAT,'X','ONCE'))
            continue
        end
        noPatoms = numel(prodAT);
        prodStoic = rxnEQTrace{i,1}(j);
        combiVect = zeros(1,noPatoms);
        for k = 1:noPatoms
            combiVect(k) = nchoosek(noPatoms,k);
        end
        combiMatrix = zeros(sum(combiVect),noPatoms);
        for k = 1:size(combiMatrix,1)
            [~,expn]=log2(max(k));
            combiMatrix(k,:) = rem(floor(k*pow2(1-max(noPatoms,expn):0)),2);
        end
        [~,index] = sort(sum(combiMatrix,2));
        combiMatrix = combiMatrix(index,:);
        
        for l = 1:size(combiMatrix,1)
            hitReactAtoms =  cell(numel(hitReact),1);
            for m = find(combiMatrix(l,:))
                for k = 1:numel(hitReact)
                    k_index = hitReact(k);
                    reactN = rxnEQTrace{i,2}{k_index};
                    reactAT = rxnEQTrace{i,3}{k_index};
                    noRatoms = numel(reactAT);
                    hitChar = regexp(reactAT, prodAT(m));
                    if ~isempty(hitChar)
                        if isempty(hitReactAtoms{k})
                            hitReactAtoms{k} = zeros(1,noRatoms);
                            hitReactAtoms{k}(hitChar) = 1;
                        else
                            hitReactAtoms{k}(hitChar) = 1;
                        end
                        break
                    end
                end
            end
            LHS = strcat([prodN, '_', char(combiMatrix(l,:)+'0')]);
            RHS = '';
            keepR = false(numel(hitReact),1);
            for m = 1:numel(hitReact)
                if ~any(hitReactAtoms{m})
                    continue
                end
                keepR(m) = true;
                reactN = rxnEQTrace{i,2}{hitReact(m)};
                if isempty(RHS)
                    RHS = strcat([reactN, '_', char(hitReactAtoms{m}+'0')]);
                else
                    RHS = strcat([RHS, ' + ', reactN, '_', char(hitReactAtoms{m}+'0')]);
                end
            end
            emuRxn = strcat([LHS,'    ',num2str(prodStoic,'%.3g'),'    v(',...
                num2str(i),')    ',RHS]);
            
            EMUrxnList{emuCC,1} = emuRxn;
            EMUrxnList{emuCC,2} = i;
            EMUrxnList{emuCC,3} = prodStoic;
            EMUrxnList{emuCC,4} = prodN;
            EMUrxnList{emuCC,5} = combiMatrix(l,:);
            EMUrxnList{emuCC,6} = rxnEQTrace{i,2}(keepR);
            EMUrxnList{emuCC,7} = hitReactAtoms(keepR);
            emuCC = emuCC + 1;
            
            %             fprintf('%s\n',emuRxn);
        end
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EMUrxnPicked = pickEMUrxn(EMUrxnList, simulatedMDVs)

EMUtracking = cell(0,2);
for i = 1:numel(simulatedMDVs)
    simMDV = strtrim(simulatedMDVs{i});
    hashPos = regexp(simMDV,'#','ONCE');
    EMUtracking{i,1} = simMDV(1:hashPos-1);
    EMUtracking{i,2} = +(simMDV(hashPos+1:end)=='1');
end

ERLsize = zeros(size(EMUrxnList,1),2);
for i = 1:size(EMUrxnList,1)
    ERLsize(i,:) = [numel(EMUrxnList{i,5}) sum(EMUrxnList{i,5})];
end

EMUrxnPicked = false(size(EMUrxnList,1),1);

while size(EMUtracking,1)>0
    EMUname = EMUtracking{1,1};
    EMUtag = EMUtracking{1,2};
    EMUlength = numel(EMUtag);
    EMUsize = sum(EMUtag);
    
    hitMet = strcmp(EMUname,EMUrxnList(:,4)) & ERLsize(:,1)==EMUlength & ...
        ERLsize(:,2)==EMUsize & ~EMUrxnPicked;
    for i = find(hitMet)'
        if ~all(EMUrxnList{i,5}==EMUtag)
            hitMet(i) = false;
        end
    end
    EMUrxnPicked(hitMet) = true;
    
    for i = find(hitMet)'
        EMUmetsHit = EMUrxnList{i,6};
        EMUtagHit = EMUrxnList{i,7};
        for j = 1:numel(EMUmetsHit)
            EMUtracking{end+1,1} = EMUmetsHit{j};
            EMUtracking{end,2} = EMUtagHit{j};
        end
    end
    
    EMUtracking(1,:) = [];
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EMUinputSubstrates = getEMUinputSubstrates (EMUrxnList)

%build input substrate list (from EMUs without reactant)
EMUreactants = cell(0,3);
for i = 1:size(EMUrxnList,1)
    for j = 1:numel(EMUrxnList{i,6})
        EMUreactants(end+1,1) = EMUrxnList{i,6}(j);
        EMUreactants(end,2) = EMUrxnList{i,7}(j);
        EMUreactants{end,3} = sum(EMUrxnList{i,7}{j});
    end
end
EMUreactants = removeEMUduplicate(EMUreactants);
hitAnyProd = false(size(EMUreactants,1),1);
for i = 1:size(EMUreactants,1)
    if any(matchEMU(EMUreactants{i,1},EMUreactants{i,2},EMUrxnList(:,[4,5])))
        hitAnyProd(i) = true;
    end
end
EMUinputSubstrates = EMUreactants(~hitAnyProd,:);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EMUout = removeEMUduplicate(EMUin)

if size(EMUin,1)==1
    EMUout = EMUin;
    return
end
cc = 1;

while cc <= size(EMUin,1)
    hitRow = strcmp(EMUin{cc,1},EMUin(:,1));
    hitRow(cc) = false;
    EMUtag = EMUin{cc,2};
    for i = find(hitRow)'
        if ~all(EMUin{i,2}==EMUtag)
            hitRow(i) = false;
        end
    end
    if any(hitRow)
        EMUin(hitRow,:) = [];
    end
    cc = cc + 1;
end
EMUout = EMUin;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hitEMU = matchEMU(EMUname, EMUtag, EMUlist)
if isempty(EMUlist)
    hitEMU = [];
    return
end

hitEMU = strcmp(EMUname,EMUlist(:,1));
for i = find(hitEMU)'
    if ~all(EMUlist{i,2}==EMUtag)
        hitEMU(i) = false;
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modelOutput = buildEMUdnsMat(EMUrxnList,EMUinputSubstrates, S, metList,matchExt)

%metList from S same as list in metListData
%get unique list of metabolites, map EMU species
%metName EMUspecies rxnIn stoicIn rxnOut stoicOut
% metListData = EMUrxnList(:,4);
% metListData = unique(metListData);
metListData = metList(~matchExt);

for i = 1:size(metListData,1)
    hitRows = strcmp(metListData{i},EMUrxnList(:,4));
    if any(hitRows)
        EMUspecies = EMUrxnList(hitRows,5);
        EMUspecies = cell2mat(EMUspecies);
        EMUspecies = unique(EMUspecies,'rows');
        
        metListData{i,2} = EMUspecies;
    else
        metListData{i,2} = [];
        EMUspecies = 0;
    end
    
    hitMet = strcmp(metListData{i},metList);
    hitInRxn = S(:,hitMet)>0;
    hitOutRxn = S(:,hitMet)<0;
    metListData{i,3} = find(hitInRxn);%rxn index
    metListData{i,4} = S(hitInRxn,hitMet);%stoic
    metListData{i,5} = find(hitOutRxn);%rxn index
    metListData{i,6} = -S(hitOutRxn,hitMet);%stoic
    metListData{i,6} = max(sum(EMUspecies,2));%molecule size
end

noIntMets = size(metListData,1);
noRxns = size(S,1);

%construct EMU and cauchy metabolites and inputs (already generated)
EMUmets = {}; %pool only single emu metabolites (unknown)
EMUmetscauchy = {}; %pool all cauchy metabolites (known)
for i = 1:size(metListData,1)
    for j = 1:size(metListData{i,2},1)
        EMUmets{end+1,1} = metListData{i,1};
        EMUmets{end,2} = metListData{i,2}(j,:);
        EMUmets{end,3} = sum(metListData{i,2}(j,:));
    end
end
for i = 1:size(EMUrxnList,1)
    if size(EMUrxnList{i,6},2)>1
        EMUmetscauchy{end+1,1} = EMUrxnList{i,6}{1};
        EMUmetscauchy{end,2} = EMUrxnList{i,7}{1};
        EMUmetscauchy{end,3} = EMUrxnList{i,6}{2};
        EMUmetscauchy{end,4} = EMUrxnList{i,7}{2};
    end
end
EMUmetscauchy = removeEMUduplicatePair(EMUmetscauchy);
for i = 1:size(EMUmetscauchy,1)
    EMUmetscauchy{i,5} = sum(EMUmetscauchy{i,2})+sum(EMUmetscauchy{i,4});
end



%construct balance matrices according to pool sizes
%use reaction map to find contributing emus
%sequence according to EMUmets
EMUsizes = unique(cell2mat(EMUmets(:,3)));%do with increment size
bigEMUmodel = cell(size(EMUsizes,1),7);

for i = 1:size(EMUsizes,1)
    hits = cell2mat(EMUmets(:,3))==EMUsizes(i);%hit single
    if isempty(EMUmetscauchy)
        hitsC = [];
    else
        hitsC = cell2mat(EMUmetscauchy(:,5))==EMUsizes(i);%hit cauchy
    end
    hitsI = cell2mat(EMUinputSubstrates(:,3))==EMUsizes(i);%hit input
    
    EMUvars = EMUmets(hits,[1,2]);
    if isempty(EMUmetscauchy)
        EMUcauchyvars = cell(0,4);
    else
        EMUcauchyvars = EMUmetscauchy(hitsC,1:4);
    end
    
    EMUinputVars = EMUinputSubstrates(hitsI,[1,2]);
    
    zeroMat = zeros(sum(hits),sum(hits) + sum(hitsI) + sum(hitsC));
    concMatMap = zeroMat; %pool size map
    fluxMatMap = cell(sum(hits),sum(hits) + sum(hitsI) + sum(hitsC));%flux map
    
    %extend conc of one met to  EMUs of same name
    for j = 1:size(EMUvars,1)
        hitMet = strcmp(EMUvars{j,1},metListData(:,1));
        concMatMap(j,j) = find(hitMet);
    end
    concMatMapList = [find(concMatMap~=0) concMatMap(concMatMap~=0)];
    %go through each EMU column and populate the flux coefficients
    %(reaction# and stoic)
    for j = 1:size(EMUvars,1)
        hitRows = matchEMU(EMUvars{j,1},EMUvars{j,2},EMUrxnList(:,[4,5]));%find EMUrxn
        for k = find(hitRows)'
            if numel(EMUrxnList{k,6})==1%match single or input
                hitReactEMU = matchEMU(EMUrxnList{k,6},EMUrxnList{k,7}{1},EMUvars);%match single
                if ~any(hitReactEMU)%match input
                    hitInputEMU = matchEMU(EMUrxnList{k,6},EMUrxnList{k,7}{1},EMUinputVars);
                    hitCol = find(hitInputEMU) + sum(hits);
                    if isempty(fluxMatMap{j,hitCol})
                        fluxMatMap{j,hitCol} = [EMUrxnList{k,2} EMUrxnList{k,3}];
                    else
                        fluxMatMap{j,hitCol} = [fluxMatMap{j,hitCol};EMUrxnList{k,2} EMUrxnList{k,3}];
                    end
                    continue
                end
                hitCol = find(hitReactEMU);
                if isempty(fluxMatMap{j,hitCol})
                    fluxMatMap{j,hitCol} = [EMUrxnList{k,2} EMUrxnList{k,3}];
                else
                    fluxMatMap{j,hitCol} = [fluxMatMap{j,hitCol};EMUrxnList{k,2} EMUrxnList{k,3}];
                end
            else%match cauchy
                EMUname1 = EMUrxnList{k,6}{1}; EMUname2 = EMUrxnList{k,6}{2};
                EMUtag1 = EMUrxnList{k,7}{1}; EMUtag2 = EMUrxnList{k,7}{2};
                hit1set1 = matchEMU(EMUname1, EMUtag1, EMUcauchyvars(:,[1 2]));
                hit1set2 = matchEMU(EMUname1, EMUtag1, EMUcauchyvars(:,[3 4]));
                hit2set1 = matchEMU(EMUname2, EMUtag2, EMUcauchyvars(:,[1 2]));
                hit2set2 = matchEMU(EMUname2, EMUtag2, EMUcauchyvars(:,[3 4]));
                hitCauchy = (hit1set1 | hit1set2) & (hit2set1 | hit2set2);
                hitCol = find(hitCauchy,1) + sum(hits) + sum(hitsI);
                
                if isempty(fluxMatMap{j,hitCol})
                    fluxMatMap{j,hitCol} = [EMUrxnList{k,2} EMUrxnList{k,3}];
                else
                    fluxMatMap{j,hitCol} = [fluxMatMap{j,hitCol};EMUrxnList{k,2} EMUrxnList{k,3}];
                end
            end
        end
    end
    
    bigEMUmodel{i,1} = EMUsizes(i);
    bigEMUmodel{i,2} = EMUvars;
    bigEMUmodel{i,3} = EMUinputVars;
    bigEMUmodel{i,4} = EMUcauchyvars;
    bigEMUmodel{i,5} = zeroMat;
    bigEMUmodel{i,6} = concMatMapList;
    bigEMUmodel{i,7} = fluxMatMap;
end

%compile unique list of flux stoic product
fluxStoicList = {};
cc = 1;
for i = 1:size(bigEMUmodel,1)
    for j = 1:size(bigEMUmodel{i,7},1)
        for k = 1:size(bigEMUmodel{i,7},2)
            if ~isempty(bigEMUmodel{i,7}{j,k})
                fluxStoicList{cc,1} = bigEMUmodel{i,7}{j,k};
                fluxStoicList{cc,2} = [i,j,k];
                cc = cc + 1;
            end
        end
    end
end
%create unique list and match indexes
fluxStoicList_unique = fluxStoicList(:,1);
delList = false(numel(fluxStoicList_unique),1);

fluxStoicList{1,3} = 1;
uniqueCounter = 2;
cr = 2;
while cr <= numel(fluxStoicList_unique)
    cc = 1;
    while cc < cr
        if delList(cc) %this is already a duplicate
            cc = cc + 1;
            continue
        end
        if size(fluxStoicList_unique{cr},1) == size(fluxStoicList_unique{cc},1)
            if size(fluxStoicList_unique{cr},1) == 1 %only 1 row
                if all(fluxStoicList_unique{cr}==fluxStoicList_unique{cc})
                    delList(cr) = true;
                    fluxStoicList{cr,3} = fluxStoicList{cc,3};
                    break
                end
            else%need to sort for 2 or more rows
                matCR = sortrows(fluxStoicList_unique{cr});
                matCC = sortrows(fluxStoicList_unique{cc});
                if all(all(matCR == matCC))
                    delList(cr) = true;
                    fluxStoicList{cr,3} = fluxStoicList{cc,3};
                    break
                end
            end
        end
        cc = cc + 1;
    end
    
    if ~delList(cr)
        fluxStoicList{cr,3} = uniqueCounter;
        uniqueCounter = uniqueCounter + 1;
    end
    cr = cr + 1;
end
fluxStoicList_unique(delList) = [];

fluxStoicT = zeros(size(fluxStoicList_unique,1),noRxns);%unique row, rxns col
for i = 1:size(fluxStoicList_unique,1)
    matNow = fluxStoicList_unique{i};
    for j = 1:size(matNow,1)
        fluxStoicT(i,matNow(j,1)) = fluxStoicT(i,matNow(j,1)) + matNow(j,2);
    end
end

for i = 1:size(bigEMUmodel,1)
    matNow = zeros(size(bigEMUmodel{i,7}));
    for j = 1:size(fluxStoicList,1)
        if fluxStoicList{j,2}(1) == i %first tag match bigEMUmodel index
            matNow(fluxStoicList{j,2}(2),fluxStoicList{j,2}(3)) = fluxStoicList{j,3};
        end
    end
    fluxMatMapList = [];
    fluxMatMapList(:,1) = find(matNow~=0);
    fluxMatMapList(:,2) = matNow(matNow~=0);
    bigEMUmodel{i,7} = fluxMatMapList;
end

modelOutput.metListData = metListData;
modelOutput.bigEMUmodel = bigEMUmodel;
modelOutput.fluxStoicT = fluxStoicT;
modelOutput.fluxStoicList_unique = fluxStoicList_unique;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EMUout = removeEMUduplicatePair(EMUin)

if size(EMUin,1)==1
    EMUout = EMUin;
    return
end
cc = 1;
while cc <= size(EMUin,1)
    name1 = EMUin{cc,1}; name2 = EMUin{cc,3};
    tag1 = EMUin{cc,2}; tag2 = EMUin{cc,4};
    
    cc_itt = cc + 1;
    while cc_itt <= size(EMUin,1)
        testname1 = EMUin{cc_itt,1}; testname2 = EMUin{cc_itt,3};
        testtag1 = EMUin{cc_itt,2}; testtag2 = EMUin{cc_itt,4};
        
        if strcmp(name1,testname1) && all(tag1==testtag1) &&...
                strcmp(name2,testname2) && all(tag2==testtag2)
            EMUin(cc_itt,:) = [];
            
        elseif strcmp(name1,testname2) && all(tag1==testtag2) &&...
                strcmp(name2,testname1) && all(tag2==testtag1)
            EMUin(cc_itt,:) = [];
        else
            cc_itt = cc_itt + 1;
        end
    end
    cc = cc + 1;
end
EMUout = EMUin;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rxnEQtext = genRxnEQ(rxnEQ)
noRxns = size(rxnEQ,1);
rxnEQtext = cell(noRxns,1);
for i = 1:noRxns
    LHS = '';
    RHS = '';
    for j = 1:numel(rxnEQ{i,1})
        if rxnEQ{i,1}(j) < 0
            if isempty(LHS)
                if abs(rxnEQ{i,1}(j)) == 1
                    LHS = rxnEQ{i,2}{j};
                else
                    LHS = strcat([num2str(abs(rxnEQ{i,1}(j))) ' ' rxnEQ{i,2}{j}]);
                end
            else
                if abs(rxnEQ{i,1}(j)) == 1
                    LHS = strcat([LHS ' + ' rxnEQ{i,2}{j}]);
                else
                    LHS = strcat([LHS ' + ' num2str(abs(rxnEQ{i,1}(j))) ' ' rxnEQ{i,2}{j}]);
                end
            end
        else
            if isempty(RHS)
                if abs(rxnEQ{i,1}(j)) == 1
                    RHS = rxnEQ{i,2}{j};
                else
                    RHS = strcat([num2str(abs(rxnEQ{i,1}(j))) ' ' rxnEQ{i,2}{j}]);
                end
            else
                if abs(rxnEQ{i,1}(j)) == 1
                    RHS = strcat([RHS ' + ' rxnEQ{i,2}{j}]);
                else
                    RHS = strcat([RHS ' + ' num2str(abs(rxnEQ{i,1}(j))) ' ' rxnEQ{i,2}{j}]);
                end
            end
        end
    end
    rxnEQtext{i} = strcat(['R' num2str(i) ' ' LHS ' = ' RHS]);
    rxnEQtext{i} = regexprep(rxnEQtext{i},'_','\\_');
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EMUinsubMID = inputSubBuilder(ofOBJ)

EMUinsubMID = cell(1,size(ofOBJ.EMUinputSubstrates,1));
modRows = false(size(ofOBJ.labelledSub,1),1);%report no match
for i = 1:size(ofOBJ.EMUinputSubstrates,1)
    ISname = ofOBJ.EMUinputSubstrates{i,1};
    emutag = ofOBJ.EMUinputSubstrates{i,2};
    hitRow = strcmp(ISname,ofOBJ.labelledSub(:,1));
    if any(hitRow) %a labelled substrate
        pos13Cenrich = ofOBJ.labelledSub{hitRow,2};
        ISmid = 1;
        for j = find(emutag==1)
            ISmid = cauchy(ISmid,cVectGen2([1-pos13Cenrich(j) pos13Cenrich(j)],1,2)')';
        end
        modRows(hitRow) = true;
    else %%not a labelled substrate
        bb6 = cVectGen2([0 1],6,7);
        ISmid = cVectGen2([1-ofOBJ.natSub13Cenrich ofOBJ.natSub13Cenrich],sum(emutag),sum(emutag)+1);
    end
    EMUinsubMID{i} = ISmid;
end
%%%report no match
if any(~modRows)
    disp('no IS match for:');
    disp(ofOBJ.labelledSub(~modRows,:));
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bigEMUmodel, cauchyTags, EMUstateStoreIS_block, EMUstate] = genEMUmodelStart(ofOBJ)
bigEMUmodel = ofOBJ.EMUmodelOutput.bigEMUmodel;
noEMUperm = size(bigEMUmodel,1);
%combine EMUcalc and EMUinSub, find numel
EMUcalcSizePerm = [];
for i = 1:noEMUperm
    EMUcombineCell = bigEMUmodel{i,2};
    EMUcalcSizePerm(i) = size(EMUcombineCell,1);
    EMUcombineCell = [EMUcombineCell; bigEMUmodel{i,3}];
    bigEMUmodel{i,8} = EMUcombineCell;
    bigEMUmodel{i,9} = size(bigEMUmodel{i,2},1);
end

%put cauchy identifiers
cauchyTags = {};
for i = 1:noEMUperm
    tagMatrix = [];
    for j = 1:size(bigEMUmodel{i,4},1)
        EMUname1 = bigEMUmodel{i,4}{j,1}; EMUname2 = bigEMUmodel{i,4}{j,3};
        EMUtag1 = bigEMUmodel{i,4}{j,2}; EMUtag2 = bigEMUmodel{i,4}{j,4};
        
        hitRow1 = cell2mat(bigEMUmodel(:,1))==sum(EMUtag1);
        hitEMUrow1 = matchEMU(EMUname1,EMUtag1,bigEMUmodel{hitRow1,8});%hit combined EMUs
        
        hitRow2 = cell2mat(bigEMUmodel(:,1))==sum(EMUtag2);
        hitEMUrow2 = matchEMU(EMUname2,EMUtag2,bigEMUmodel{hitRow2,8});
        
        tagMatrix(end+1,:) = [find(hitRow1) find(hitEMUrow1) find(hitRow2) find(hitEMUrow2)];
    end
    cauchyTags{end+1,1} = tagMatrix;
end
%%%EMU inputsubstrates and initial
for i = 1:noEMUperm
    for j = 1:size(bigEMUmodel{i,3},1)
        hit = matchEMU(bigEMUmodel{i,3}{j,1},bigEMUmodel{i,3}{j,2}, ofOBJ.EMUinputSubstrates);
        bigEMUmodel{i,3}{j,3} = find(hit);
    end
end
%             EMUinsubMID = EMUModel_substrate_mod_whole();%modify this with new model
EMUstateStoreIS = cell(ofOBJ.noSteps,size(ofOBJ.EMUinputSubstrates,1));
for i = 1:size(ofOBJ.EMUinputSubstrates,1)
    for j = 1:ofOBJ.noSteps
        EMUstateStoreIS{j,i} = ofOBJ.EMUinsubMID{i};
    end
end
%reparse input substrates into blocks
EMUstateStoreIS_block = cell(ofOBJ.noSteps,noEMUperm);
for i = 1:noEMUperm
    if ~isempty(bigEMUmodel{i,3})
        blockRow = cell2mat(bigEMUmodel{i,3}(:,3));
        for j = 1:ofOBJ.noSteps
            EMUstateStoreIS_block{j,i} = cell2mat(EMUstateStoreIS(j,blockRow))';
        end
    end
end

%EMUtoBeCalculated & EMUinputSubtrate, EMUcauchy
EMUstate = cell(noEMUperm,2);

for i = 1:noEMUperm
    EMUsize = bigEMUmodel{i,1};
    EMUvect = cVectGen2([1-ofOBJ.natEndo13Cenrich, ofOBJ.natEndo13Cenrich],EMUsize,EMUsize+1)';
    matNowCalc = EMUvect(ones(size(bigEMUmodel{i,2},1),1),:);
    matNowCauchy = EMUvect(ones(size(bigEMUmodel{i,4},1),1),:);
    EMUstate{i,1} = [matNowCalc; EMUstateStoreIS_block{1,i}];
    EMUstate{i,2} = matNowCauchy;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [noSteps, simTime, tSampleIndex, noParaPerFlux, knotSeq, knotSeqInt, Nout, Nout_int] = genBSplineMat(ofOBJ)
%generate BSpline basis matrix from knot and step
tSimScale = [];
tSampleIndex = [];
stepBox = ofOBJ.stepBTWsample;
tSample = ofOBJ.sampleTime;
tScale = tSample/max(tSample);
for i = 1:numel(tScale)
    if i == 1
        tStep = (tScale(i) - 0)/stepBox(i);
        tSimScale = [tSimScale 0:tStep:tScale(i)];
        tSampleIndex(i) = numel(tSimScale);
        continue
    end
    tStep = (tScale(i) - tScale(i-1))/stepBox(i);
    tSimScale = [tSimScale tScale(i-1)+tStep:tStep:tScale(i)];
    tSampleIndex(i) = numel(tSimScale);
end
tSimNoScale = tSimScale*max(tSample);
deltaT = [0 diff(tSimScale)]*max(tSample);
noSteps = numel(tSimScale);
simTime = tSimNoScale;

orderS = ofOBJ.orderS;
noKnots = ofOBJ.noIntKnots+2*orderS;
noParaPerFlux = noKnots - orderS;
xKnot = ofOBJ.intKntPos;
knotSeq = [zeros(1,orderS) xKnot ones(1,orderS)];
Nout = bSplineMat(knotSeq,tSimScale,orderS);
knotSeqInt = [knotSeq 1];%%didn't add a LHS zero (then don't need [0 CP])
Nout_int = bSplineMat(knotSeqInt,tSimScale,orderS+1);
tIntMat = zeros(noParaPerFlux);
for i = 1:noParaPerFlux
    tIntMat(i,i:end) = knotSeqInt(i+orderS)-knotSeqInt(i);
end
Nout_int = tIntMat*tSample(end)/orderS*Nout_int;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dataOut = readMetDatFile(fileName)
fid = fopen(fileName,'r');
c = textscan(fid,'%s','delimiter','\n');
fclose(fid);
c = c{1};
%%%%identify data blocks, remove initial blanks
while isempty(c{1})
    c(1) = [];
end

contiBlock = [];%start,stop inclusive
rc = 1;
startRC = rc;
scanStop = true;
while rc<=size(c,1)
    if scanStop && isempty(c{rc})
        contiBlock(end+1,:) = [startRC rc-1];
        scanStop = false;
    elseif ~scanStop && ~isempty(c{rc})
        startRC = rc;
        scanStop = true;
    end
    rc = rc + 1;
end
dataOut = cell(size(contiBlock,1),3);
for i = 1:size(contiBlock,1)
    dataOut{i,1} = strtrim(c{contiBlock(i,1)});
    dataBlock = c(contiBlock(i,1)+1:contiBlock(i,2));
    dataBlockVal = [];
    for j = 1:size(dataBlock,1)
        tabPos = regexp(dataBlock{j},'\t');
        tabPos = [0 tabPos];
        cc = 1;
        for k = 1:numel(tabPos)-1
            dataBlockVal(j,cc) = str2double(dataBlock{j}(tabPos(k)+1:tabPos(k+1)-1));
            cc = cc + 1;
        end
        dataBlockVal(j,cc) = str2double(dataBlock{j}(tabPos(k)+1:end));
    end
    dataBlockVal(isnan(dataBlockVal)) = 0;
    dataOut{i,2} = dataBlockVal(:,1:2:size(dataBlockVal,2));
    dataOut{i,3} = dataBlockVal(:,2:2:size(dataBlockVal,2));
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dataOut = readIonFormFile(fileName)
fid = fopen(fileName,'r');
c = textscan(fid,'%s','delimiter','\n');
fclose(fid);
c = c{1};
%%%%identify data blocks, remove initial blanks
while isempty(c{1})
    c(1) = [];
end
dataOut = cell(numel(c),3);
for i = 1:numel(c)
    tabPos = regexp(c{i},'\t');
    dataOut{i,1} = strtrim(c{i}(1:tabPos(1)-1));
    emuTag = c{i}(tabPos(1)+1:tabPos(2)-1);
    emuTagMat = [];
    for j = 1:numel(emuTag)
        emuTagMat(j) = str2double(emuTag(j));
    end
    dataOut{i,2} = emuTagMat;
    dataOut{i,3} = strtrim(c{i}(tabPos(2)+1:end));
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dataMet = calcErrorByMC(dataIn,noItt)
for i = 1:size(dataIn,1)
    fracAb = dataIn{i,2};
    fractAbSE = dataIn{i,3};
    sizeAB = size(fracAb);
    %%%need to populate zero values before calculating fractions%%%
    fracAb(fracAb<0) = 0;
    totalMat = sum(fracAb,1);
    totalMatSE = zeros(1,sizeAB(2));
    fractMat = fracAb./totalMat(ones(sizeAB(1),1),:);
    %%%%%%%%
    fractMatSE = zeros(size(fracAb));
    for j = 1:sizeAB(2)
        abUse = fracAb(:,j);
        errUse = fractAbSE(:,j);
        allMat = zeros(sizeAB(1),noItt);
        for k = 1:noItt
            allMat(:,k) = abUse + errUse.*randn(sizeAB(1),1);
        end
        allMat(allMat<0) = 0;
        allMat_sum = sum(allMat,1);
        allMat_fract = allMat./allMat_sum(ones(sizeAB(1),1),:);
        
        totalMatSE(j) = std(allMat_sum);
        hitZero = allMat_sum ==0;
        fractMatSE(:,j) = std(allMat_fract(:,~hitZero)')';
    end
    if any(totalMat==0)
        hitCol = totalMat==0;
        fractMat(:,hitCol) = 0;
        fractMatSE(:,hitCol) = 0;
    end
    dataIn{i,4} = totalMat;%total witout empty value fill
    dataIn{i,5} = totalMatSE;%error estimate with empty value fill
    dataIn{i,6} = fractMat;%fraction with empty value fill
    dataIn{i,7} = fractMatSE;%error estimate with empty value fill
end
dataMet = dataIn;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dataMet = expandDataMet(ofOBJ)
dataMet = ofOBJ.dataMet;
expMIDs = ofOBJ.ionForm;
noSamples = numel(ofOBJ.sampleTime);
cc = 0;
for i = 1:size(dataMet,1)
    dataMetName = dataMet{i,1};
    
    fractKeep = dataMet{i,7}~=0;%this will eliminate fract of 1
    %remove 1 (most bottom) dependent data from normalization of a column
    for j = 1:noSamples
        hitRow = find(fractKeep(:,j),1,'last');
        if ~isempty(hitRow)
            fractKeep(hitRow,j) = false;
        end
    end
    fractKeep = find(fractKeep);
    hitIntMet = strcmp(dataMet{i,1},ofOBJ.metListInt);
    
    formulaRow = find(strcmp(dataMetName,expMIDs(:,1)));
    metEMUsize =  sum(expMIDs{formulaRow,2});
    hitRow = cell2mat(ofOBJ.bigEMUmodel(:,1))==metEMUsize;
    hitEMU = matchEMU(expMIDs{formulaRow,1},expMIDs{formulaRow,2},ofOBJ.bigEMUmodel{hitRow,2});
    CM = corrMatGen2([1:metEMUsize+1],metEMUsize+1,expMIDs{formulaRow,3});%this one is not normalized, for abs calc
    
    unlabelledVect =  ofOBJ.EMUstate{hitRow,1}(hitEMU,:)';
    
    dataMet{i,8} = find(hitIntMet);
    dataMet{i,9} = [fractKeep [1:numel(fractKeep)]'+cc];
    dataMet{i,10} = [find(hitRow) find(hitEMU)];
    dataMet{i,11} = CM;%
    dataMet{i,12} = metEMUsize+1;%EMU size
    dataMet{i,13} = unlabelledVect(:,ones(1,noSamples));
    dataMet{i,14} = dataMet{i,6}==0;%discard away zero points
    
    cc = cc+numel(fractKeep);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opCon, opInput] = genOptimisationProblem(ofOBJ)
opCon = []; opInput =[];
noExpData = size(ofOBJ.dataMet,1);
noSamples = numel(ofOBJ.sampleTime);
%vectorize exp data
dataMetConc_EXP = cell2mat(ofOBJ.dataMet(:,4));
dataMetSE_EXP = cell2mat(ofOBJ.dataMet(:,5));
dataMetConc_EXPvect = dataMetConc_EXP(:);
dataMetSE_EXPvect = dataMetSE_EXP(:);
metConcProfile_SIM = zeros(noExpData,noSamples);

midTotLength = ofOBJ.dataMet{end,9}(end,2);
dataMetMID_EXPvect = zeros(midTotLength,1);
dataMetMID_SEvect = zeros(midTotLength,1);
dataMetMID_SIMvect = dataMetMID_EXPvect;%clone
for i = 1:size(ofOBJ.dataMet,1)
    dataMetMID_EXPvect(ofOBJ.dataMet{i,9}(:,2)) = ofOBJ.dataMet{i,6}(ofOBJ.dataMet{i,9}(:,1));
    dataMetMID_SEvect(ofOBJ.dataMet{i,9}(:,2)) = ofOBJ.dataMet{i,7}(ofOBJ.dataMet{i,9}(:,1));
end
dataMetMID_SEvect(dataMetMID_SEvect<ofOBJ.midMinError) = ofOBJ.midMinError;%set minimum error to 1%

f_base = ofOBJ.fluxScale(:,ones(ofOBJ.noParaPerFlux,1));
f_diff = ofOBJ.fluxScale(:,2)-ofOBJ.fluxScale(:,1);
f_diff = f_diff(:,ones(ofOBJ.noParaPerFlux,1));

c_base = ofOBJ.concScale(:,1);%this is for the initial concentration
c_diff = ofOBJ.concScale(:,2)-ofOBJ.concScale(:,1);
cstag = ofOBJ.concScale(cell2mat(ofOBJ.dataMet(:,8)),2);%take maximum

noX = ofOBJ.noParaPerFlux*ofOBJ.noReactions+ofOBJ.noIntMets+noExpData;
lb = zeros(noX,1);%organized in rxn blocks, then met
ub = ones(noX,1);

%construct flux profile from bspline basis vector
CPmap = zeros(ofOBJ.noParaPerFlux,ofOBJ.noReactions);
CPmap(:) = 1:ofOBJ.noParaPerFlux*ofOBJ.noReactions;
CPmap = CPmap';

%construct initial metabolite
concMap = [1:ofOBJ.noIntMets] + ofOBJ.noParaPerFlux*ofOBJ.noReactions;%met conc comes after fluxes
%construct metabolite active fraction
stagMap = [1:noExpData] + ofOBJ.noParaPerFlux*ofOBJ.noReactions + ofOBJ.noIntMets;

f_diff_vect = f_diff';
diff_vect = [f_diff_vect(:);c_diff];
diff_vect(noX) = 0;
f_base_vect = f_base';
base_vect = [f_base_vect(:);c_base];
base_vect(noX) = 0;

%%%n mat incorporate stoichiometry%%
bigNmat = zeros(ofOBJ.noIntMets*ofOBJ.noSteps,noX);
tIndex = 1:ofOBJ.noIntMets;
for tSlice = 1:ofOBJ.noSteps
    bigN = zeros(ofOBJ.noIntMets,noX);
    for i = 1:ofOBJ.noIntMets
        NintSlice = ofOBJ.Nout_int(:,tSlice)';
        for j = 1:ofOBJ.noReactions
            Ncols = CPmap(j,:);
            bigN(i,Ncols) = NintSlice*ofOBJ.Sint(i,j)+bigN(i,Ncols);
        end
        bigN(i,concMap(i)) = 1;
    end
    bigNmat((tSlice-1)*ofOBJ.noIntMets+tIndex,:) = bigN;
end
Acon = sparse(-bigNmat*diag(diff_vect));
Bcon = bigNmat*base_vect-ofOBJ.concBound(1)*ones(ofOBJ.noSteps*ofOBJ.noIntMets,1);

opCon.Acon = Acon;
opCon.Bcon = Bcon;
opCon.lb = lb;
opCon.ub = ub;

opInput.CPmap = CPmap;
opInput.concMap = concMap;
opInput.f_base = f_base;
opInput.f_diff = f_diff;
opInput.c_base = c_base;
opInput.c_diff = c_diff;
opInput.noT = ofOBJ.noSteps;
opInput.Nout_int = ofOBJ.Nout_int;
opInput.noEMUperm = size(ofOBJ.bigEMUmodel,1);
opInput.cauchyTags = ofOBJ.cauchyTags;
opInput.Sint = ofOBJ.Sint;
opInput.deltaT = [0 diff(ofOBJ.simTime)];
opInput.tSampleIndex = ofOBJ.tSampleIndex;
opInput.noTsample = numel(ofOBJ.tSampleIndex);
opInput.dataMet = ofOBJ.dataMet;
opInput.fluxStoicT = ofOBJ.EMUmodelOutput.fluxStoicT;
opInput.Nout = ofOBJ.Nout;
opInput.noExpData = noExpData;
opInput.dataMetConc_EXPvect = dataMetConc_EXPvect;
opInput.dataMetSE_EXPvect = dataMetSE_EXPvect;
opInput.dataMetMID_SIMvect = dataMetMID_SIMvect;
opInput.dataMetMID_EXPvect = dataMetMID_EXPvect;
opInput.dataMetMID_SEvect = dataMetMID_SEvect;
opInput.EMUstateStoreIS_block = ofOBJ.EMUstateStoreIS_block;
opInput.EMUstateStore = cell(ofOBJ.noSteps,1);
opInput.EMUstate = ofOBJ.EMUstate;
opInput.stagMap = stagMap;
opInput.cstag = cstag;
opInput.metConcProfile_SIM = metConcProfile_SIM;
opInput.A_cell = ofOBJ.bigEMUmodel(:,5);
opInput.Cmap_cell = ofOBJ.bigEMUmodel(:,6);
opInput.Vmap_cell = ofOBJ.bigEMUmodel(:,7);
opInput.EMUsize = cell2mat(ofOBJ.bigEMUmodel(:,1))+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function genConcScale(ofOBJ)
for i = 1:size(ofOBJ.dataMet)
    minimumConc = min(ofOBJ.dataMet{i,4} - 4*ofOBJ.dataMet{i,5});
    if minimumConc < 0
        minimumConc = ofOBJ.concBound(1);
    end
    maximumConc = max(ofOBJ.dataMet{i,4} + 4*ofOBJ.dataMet{i,5});
    hitIntMet = strcmp(ofOBJ.dataMet{i,1},ofOBJ.metListInt);
    ofOBJ.concScale(hitIntMet,:) = [minimumConc maximumConc];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%