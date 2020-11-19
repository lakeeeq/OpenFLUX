function drawFig3(source,callbackdata,mid_outParsed,fMID_diff, OF,flux)
%for steady state
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

pp1 = subplot(2,1,1);

ylabel('MID');
if isempty(fMID_diff)
    h1 = bar([mid_outParsed{toDoF,2}';mid_outParsed{toDoF,2}'*0],'stacked');
    pp1.XTickLabel = {'sim',''};
else
    h1 = bar([mid_outParsed{toDoF,2}';mid_outParsed{toDoF,3}'],'stacked');
    title(strcat([metName ' ' num2str(mid_outParsed{toDoF,5})]));
    pp1.XTickLabel = {'sim','exp'};
end

if ~isempty(fMID_diff)
    hold on
    dataMIDAve = mid_outParsed{toDoF,3};
    dataMIDSE = mid_outParsed{toDoF,4};
    for j = 1:midSize
        lowB = sum(dataMIDAve(1:j)) - 2*dataMIDSE(j);
        highB = sum(dataMIDAve(1:j)) + 2*dataMIDSE(j);
        plot([2 2],[lowB highB], 'k-','LineWidth',10);
    end
    hold off
end
legend(massFract);

pp2 = subplot(2,1,2);
h2 = bar(flux);
pp2.YScale = 'log';
ylabel('flux');
xlabel('reactions');

drawFigFH = @(sourceX,callbackdataX)drawFig3(sourceX,callbackdataX,...
    mid_outParsed,fMID_diff, OF,flux);

PBn = uicontrol('Style', 'pushbutton','String','next','Position',[80 20 60 20],...
    'Callback',drawFigFH,'Value',toDoF,'Max',toDoF);
PBb = uicontrol('Style', 'pushbutton','String','back','Position',[20 20 60 20],...
    'Callback',drawFigFH,'Value',toDoF,'Max',toDoF);
