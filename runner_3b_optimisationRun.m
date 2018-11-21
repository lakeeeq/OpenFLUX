function runner_3b_optimisationRun(fileToLoad)
opSaveFolder = 'OPinstances/';
disp(['loading file ' fileToLoad]);
load(strcat(opSaveFolder,fileToLoad));
disp([opSave.saveFileName ' loaded']);
opOptions = optimoptions('fmincon', 'Display','iter','MaxFunEvals',40000);

if ~opSave.isODE
    [Acon, Bcon] = OpenFLUX.buildSBRcon(opSave.AconParas);
end

runFeas = false;
if isempty(opSave.xFeas)%get first feasible soln
    runFeas = true;
    x0 = opSave.x0;
    xGuess = x0;
elseif opSave.solnIndex == 0%%there is a feasible soln, but not optimised
    xStart = opSave.xFeas;
elseif ~opSave.isODE && ~all(opSave.stepBTWsample==opSave.xFitSeries(opSave.solnIndex).stepBTWsample) %change in resolution
    runFeas = true;
    x0 = opSave.xFitSeries(opSave.solnIndex).xFinish;%not too far from last op soln
    xGuess = x0;
elseif opSave.isODE && ~isempty(opSave.xFitSeries(opSave.solnIndex).stepBTWsample) %change from SBR to ODE
    runFeas = true;
    x0 = opSave.xFitSeries(opSave.solnIndex).xFinish;%not too far from last op soln
    xGuess = x0;
else %%get latest feasible soln
    xStart = opSave.xFitSeries(opSave.solnIndex).xFinish;
end

if runFeas
    disp('finding/checking for feasible soln...')
    if opSave.isODE
        while 1
            xFeas = fmincon(@(x)minDistX0(x,x0),xGuess,[],[],[],[],opSave.lb,opSave.ub,opSave.conFxn,opOptions);
            if all(opSave.conFxn(xFeas)<=0)
                break
            else
                xGuess = xFeas;
            end
        end
    else
        while 1
            xFeas = fmincon(@(x)minDistX0(x,x0),xGuess,Acon,Bcon,[],[],opSave.lb,opSave.ub,[],opOptions);
            if all(Acon*xFeas<=Bcon)
                break
            else
                xGuess = xFeas;
            end
        end
    end
    xStart = xFeas;
    if isempty(opSave.xFeas)
        opSave.xFeas = xFeas;
        save(strcat(opSaveFolder, opSave.saveFileName), 'opSave');
    end
end


tStart = tic;
disp('running optimisation, once...')
if opSave.isODE
    [xFinish,fval,exitflag]= fmincon(opSave.fitFxn,xStart,[],[],[],[],opSave.lb,opSave.ub,opSave.conFxn,opOptions);
%     xFinish = xStart;fval = 1e11;exitflag = -1;
else   
    [xFinish,fval,exitflag]= fmincon(opSave.fitFxn,xStart,Acon,Bcon,[],[],opSave.lb,opSave.ub,[],opOptions);
%     xFinish = xStart;fval = 1e11;exitflag = -1; 
end
tElapse = toc(tStart);

opSave.solnIndex = opSave.solnIndex+1;
opSave.xFitSeries(opSave.solnIndex).saveTime = datetime;
opSave.xFitSeries(opSave.solnIndex).xStart = xStart;
opSave.xFitSeries(opSave.solnIndex).xFinish = xFinish;
opSave.xFitSeries(opSave.solnIndex).tElapse = tElapse;
opSave.xFitSeries(opSave.solnIndex).fval = fval;
opSave.xFitSeries(opSave.solnIndex).exitflag = exitflag;
if ~opSave.isODE
    opSave.xFitSeries(opSave.solnIndex).stepBTWsample = opSave.stepBTWsample;
end
save(strcat(opSaveFolder, opSave.saveFileName), 'opSave');
disp('optimisation saved...')