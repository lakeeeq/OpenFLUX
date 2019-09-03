function runner_3b_optimisationRun(opSaveFolder,fileToLoad)
cauchy(1,[1 0]);
disp(['loading file ' fileToLoad]);
load(strcat(opSaveFolder,fileToLoad));
disp([opSave.saveFileName ' loaded']);
opOptions = optimoptions('fmincon', 'Display','iter','MaxFunEvals',40000);

if ~opSave.isODE && opSave.isDynamic
    [Acon, Bcon] = OpenFLUX.buildSBRcon(opSave.AconParas);%build locally because big matrices
end

runFeas = false;
if isempty(opSave.xFeas)%get first feasible soln
    runFeas = true;
    x0 = opSave.x0;
    xGuess = x0;
elseif opSave.solnIndex == 0%%there is a feasible soln, but not optimised
    xStart = opSave.xFeas;
elseif ~opSave.isODE && opSave.isDynamic && ~all(opSave.stepBTWsample==opSave.xFitSeries(opSave.solnIndex).stepBTWsample) %change in SBR step size
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
        if opSave.bumpUpXFeasIniConc
            xFeas = OF.simSolnODE_stepConc(xFeas);
        end
    elseif opSave.isDynamic
        while 1
            xFeas = fmincon(@(x)minDistX0(x,x0),xGuess,Acon,Bcon,[],[],opSave.lb,opSave.ub,[],opOptions);
            if all(Acon*xFeas<=Bcon)
                break
            else
                xGuess = xFeas;
            end
        end
    else
        while 1
            xFeas = fmincon(@(x)minDistX0(x,x0),xGuess,[],[],[],[],opSave.lb,opSave.ub,opSave.conFxn,opOptions);
            if all(opSave.conFxn(xFeas)<=0)
                break
            else
                xGuess = xFeas;
            end
        end
    end
    xStart = xFeas;
    if isempty(opSave.xFeas)
        opSave.xFeas = xFeas;
        save(strcat(opSaveFolder, opSave.saveFileName), 'opSave','OF');
    end
end


tStart = tic;
disp('running optimisation, once...')
if opSave.isODE%ODE mode
    [xFinish,fval,exitflag]= fmincon(opSave.fitFxn,xStart,[],[],[],[],opSave.lb,opSave.ub,opSave.conFxn,opOptions);
    %     xFinish = xStart;fval = 1e11;exitflag = -1;
elseif opSave.isDynamic%SBR mode
    [xFinish,fval,exitflag]= fmincon(opSave.fitFxn,xStart,Acon,Bcon,[],[],opSave.lb,opSave.ub,[],opOptions);
    %     xFinish = xStart;fval = 1e11;exitflag = -1;
else%SS mode
    [xFinish,fval,exitflag]= fmincon(opSave.fitFxn,xStart,[],[],[],[],opSave.lb,opSave.ub,opSave.conFxn,opOptions);
end
tElapse = toc(tStart);

opSave.solnIndex = opSave.solnIndex+1;
opSave.xFitSeries(opSave.solnIndex).saveTime = datetime;
opSave.xFitSeries(opSave.solnIndex).xStart = xStart;
opSave.xFitSeries(opSave.solnIndex).xFinish = xFinish;
opSave.xFitSeries(opSave.solnIndex).tElapse = tElapse;
opSave.xFitSeries(opSave.solnIndex).fval = fval;
opSave.xFitSeries(opSave.solnIndex).exitflag = exitflag;
if ~opSave.isODE && opSave.isDynamic
    opSave.xFitSeries(opSave.solnIndex).stepBTWsample = opSave.stepBTWsample;
end
save(strcat(opSaveFolder, opSave.saveFileName), 'opSave','OF');
disp('optimisation saved...');
disp(strcat(opSaveFolder, opSave.saveFileName));