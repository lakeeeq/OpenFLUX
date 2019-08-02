function drawFig1(source,callbackdata,mid_outParsed,tSimNoScale,tSample)

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
if toDoF > size(mid_outParsed,1)
    toDoF = size(mid_outParsed,1);
end
if toDoF < 1
    toDoF = 1;
end

hsvMap = colormap('hsv');
hsvMap = hsvMap(ceil(1:64/8:64),:);
hsvMapflux = colormap('hsv');
hsvMapflux = hsvMapflux(ceil(1:64/15:64),:);
midSize = size(mid_outParsed{toDoF,2},1);
massFract = {};
for j = 1:midSize
    massFract{j} = strcat(['m' num2str(j-1)]);
end
metName = mid_outParsed{toDoF,1};

subplot(5,1,1);
h0 = plot(tSimNoScale, sum(mid_outParsed{toDoF,2},1));
title(metName);
ylabel('met conc.');
subplot(5,1,[2 3]);
h1 = plot(tSimNoScale, mid_outParsed{toDoF,2});
for j = 1:numel(h1)
    set(h1(j),'Color',hsvMap(j,:));
end
if ~isempty(mid_outParsed{toDoF,7})
    hold on
    dataMIDAve = mid_outParsed{toDoF,7};
    dataMIDSE = mid_outParsed{toDoF,8};
    plotMat = dataMIDSE~=0;
    for j = 1:midSize
        plot(tSample(plotMat(j,:)), dataMIDAve(j,plotMat(j,:)),'x','Color',hsvMap(j,:));
        for k = find(plotMat(j,:))
            lowB = dataMIDAve(j,k) - 2*dataMIDSE(j,k);
            highB = dataMIDAve(j,k) + 2*dataMIDSE(j,k);
            plot([tSample(k) tSample(k)],[lowB highB], '-','Color',hsvMap(j,:));
        end

    end
    hold off
end
ylabel('isotop. conc.');
legend(massFract);

rxnIn = mid_outParsed{toDoF,3};
noRxnIn = size(rxnIn,1);
rxnOut = mid_outParsed{toDoF,4};
vSub = [rxnIn;rxnOut];
vLegend = [mid_outParsed{toDoF,5};mid_outParsed{toDoF,6}];

subplot(5,1,[4 5]);
h2 = plot(tSimNoScale, vSub);
for j = 1:noRxnIn
    set(h2(j),'Color',hsvMapflux(j,:));
end
for j = 1+noRxnIn:size(vSub,1)
    set(h2(j),'Color',hsvMapflux(j,:),'LineStyle','--');
end
ylabel('fluxes');
xlabel('time(min)');
legend(vLegend);

drawFigFH = @(sourceX,callbackdataX)drawFig1(sourceX,callbackdataX,...
    mid_outParsed,tSimNoScale,tSample);

PBn = uicontrol('Style', 'pushbutton','String','next','Position',[80 20 60 20],...
    'Callback',drawFigFH,'Value',toDoF,'Max',toDoF);
PBb = uicontrol('Style', 'pushbutton','String','back','Position',[20 20 60 20],...
    'Callback',drawFigFH,'Value',toDoF,'Max',toDoF);
