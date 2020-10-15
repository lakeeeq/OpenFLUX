%%%create OpenFLUX object and read txt model%%%
%%%specify input substrate%%%
%%%specify time of samples%%%
%%%specifcy ion formula file%%% 

OFspec = OpenFLUX.OFobjSpecification(mfilename,OFspecFileName);

OF = OpenFLUX;
% return
%%%%initial configuration%%%%
OF.modelCondition = OFspec.modelCondition;
OF.modelFileName = OFspec.modelFileName;
OF.natEndo13Cenrich  = OFspec.natEndo13Cenrich;
OF.natSub13Cenrich = OFspec.natSub13Cenrich;
OF.fluxBound = OFspec.fluxBound;
OF.isDynamic = OFspec.isDynamic;
if OF.isDynamic
    OF.concBound = OFspec.concBound;
else
    OF.par = OFspec.par;
    OF.ufISA = OFspec.unlabelledFraction;
end
        
OF.buildModel;
OF.labelledSub = OFspec.labelledSub;

if OF.isDynamic
    OF.intKntPos = OFspec.intKntPos;
    OF.orderS = OFspec.orderS;
    OF.sampleTime = OFspec.sampleTime;
    OF.odeSimTime = OFspec.odeSimTime;
    OF.isODEsolver = OFspec.isODEsolver;
    if ~OF.isODEsolver
        OF.stepBTWsample = OFspec.stepBTWsample;
    end
end

save(strcat(OFspec.inputDirectory,filesep,OFspec.modelObjSaveName),'OF');

emuStateMap = OF.getPrivProp('emuMapping');
fprintf('list of EMUs and positions:\n');
fprintf('metName(col2)\t[EMU](col3)\t[EMUindexMap](col4)\n');
for i = 1:size(OF.emuList,1)
    fprintf('%1.0f\t%s\t[%s]\t[%1.0f,%1.0f]\n',i,OF.emuList{i,1},regexprep(num2str(OF.emuList{i,2}),'  ',','),emustateMap(i,1),emustateMap(i,2));
end