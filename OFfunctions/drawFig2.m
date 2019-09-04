function drawFig1(source,callbackdata,mid_outParsed_ref,compiledResults,simTime_ref,sampleTime_ref)

if isnumeric(source) && ~isempty(source)
    toDoF = source;
elseif ~isempty(source)
    toDoF = source.Value;
    if strcmp(source.String,'next')
        toDoF = toDoF + 1;
    else
        toDoF = toDoF - 1;
    end
else
    toDoF = 1;
end
if toDoF > size(mid_outParsed_ref,1)
    toDoF = size(mid_outParsed_ref,1);
end
if toDoF < 1
    toDoF = 1;
end

hsvMap = colormap('hsv');
hsvMap = hsvMap(ceil(1:64/8:64),:);
hsvMapflux = colormap('hsv');
hsvMapflux = hsvMapflux(ceil(1:64/15:64),:);
midSize = size(mid_outParsed_ref{toDoF,2},1);
massFract = {};
for j = 1:midSize
    massFract{j} = strcat(['m' num2str(j-1)]);
end
metName = mid_outParsed_ref{toDoF,1};
f2 = figure(2);
clf(f2);
f2.Position = [809 263 696 666];
subplot(5,1,1);
hold on
for i = 1:numel(compiledResults)
    h0 = plot(simTime_ref, sum(compiledResults(i).mid_outParsed{toDoF,2},1),'-k');
end
hold off
title(metName);
ylabel('met conc.');
subplot(5,1,[2 3]);
hold on
for i = 1:numel(compiledResults)
    h1 = plot(simTime_ref, compiledResults(i).mid_outParsed{toDoF,1});
    for j = 1:numel(h1)
        set(h1(j),'Color',hsvMap(j,:));
    end
end
hold off

if ~isempty(mid_outParsed_ref{toDoF,7})
    hold on
    dataMIDAve = mid_outParsed_ref{toDoF,7};
    dataMIDSE = mid_outParsed_ref{toDoF,8};
    plotMat = dataMIDSE~=0;
    for j = 1:midSize
        plot(sampleTime_ref(plotMat(j,:)), dataMIDAve(j,plotMat(j,:)),'x','Color',hsvMap(j,:));
        for k = find(plotMat(j,:))
            lowB = dataMIDAve(j,k) - 2*dataMIDSE(j,k);
            highB = dataMIDAve(j,k) + 2*dataMIDSE(j,k);
            plot([sampleTime_ref(k) sampleTime_ref(k)],[lowB highB], '-','Color',hsvMap(j,:));
        end

    end
    hold off
end
ylabel('isotop. conc.');
legend(massFract);

subplot(5,1,[4 5]);
vLegend = [mid_outParsed_ref{toDoF,5};mid_outParsed_ref{toDoF,6}];
hold on
for i = 1:numel(compiledResults)
    rxnIn = compiledResults(i).mid_outParsed{toDoF,2};
    noRxnIn = size(rxnIn,1);
    rxnOut = compiledResults(i).mid_outParsed{toDoF,3};
    vSub = [rxnIn;rxnOut];
    
    h2 = plot(simTime_ref, vSub);
    for j = 1:noRxnIn
        set(h2(j),'Color',hsvMapflux(j,:));
    end
    for j = 1+noRxnIn:size(vSub,1)
        set(h2(j),'Color',hsvMapflux(j,:),'LineStyle','--');
    end
end
hold off
ylabel('fluxes');
xlabel('time (min)');
legend(vLegend);

drawFigFH = @(sourceX,callbackdataX)drawFig2(sourceX,callbackdataX,...
    mid_outParsed_ref,compiledResults,simTime_ref,sampleTime_ref);

PBn = uicontrol('Style', 'pushbutton','String','next','Position',[80 20 60 20],...
    'Callback',drawFigFH,'Value',toDoF,'Max',toDoF);
PBb = uicontrol('Style', 'pushbutton','String','back','Position',[20 20 60 20],...
    'Callback',drawFigFH,'Value',toDoF,'Max',toDoF);
