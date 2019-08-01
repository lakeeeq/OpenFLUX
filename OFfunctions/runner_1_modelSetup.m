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
OF.ionFormFileName = OFspec.ionFormFileName;
OF.natEndo13Cenrich  = OFspec.natEndo13Cenrich;
OF.natSub13Cenrich = OFspec.natSub13Cenrich;
OF.concBound = OFspec.concBound;
OF.fluxBound = OFspec.fluxBound;
        
OF.buildModel;

OF.intKntPos = OFspec.intKntPos;
OF.orderS = OFspec.orderS;
OF.labelledSub = OFspec.labelledSub;
OF.sampleTime = OFspec.sampleTime;
OF.stepBTWsample = OFspec.stepBTWsample;
OF.odeSimTime = OFspec.odeSimTime;
OF.isODEsolver = OFspec.isODEsolver;
save(OFspec.modelObjSaveName,'OF');