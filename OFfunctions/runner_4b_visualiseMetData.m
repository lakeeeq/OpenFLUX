%%%plot exp data contained in OpenFLUX instances%%%

OFspec = OpenFLUX.OFobjSpecification(mfilename,OFspecFileName);
load(strcat(OFspec.loadFolder,OFspec.fileName));

if ~OF.isOptimisation
    disp('only optimisation instances contain exp data. terminate.');
    return
end
disp(['modelTag: ' OF.modelCondition]);
disp(['dataTag: ' OF.metDataFileName]);

drawFig0([],[],OF.dataMet,OF.sampleTime);