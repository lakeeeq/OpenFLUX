function drawFig0(source,callbackdata,dataMet,t)

if ~isempty(source)
    toDoF = source.Value;
    if strcmp(source.String,'next')
        toDoF = toDoF + 1;
    else
        toDoF = toDoF - 1;
    end
else
    toDoF = 1;
end
if toDoF > size(dataMet,1)
    toDoF = size(dataMet,1);
end
if toDoF < 1
    toDoF = 1;
end

hsvMap = colormap('hsv');
hsvMap = hsvMap(ceil(1:64/10:64),:);
midSize = size(dataMet{toDoF,2},1);
massFract = cell(midSize,1);
for j = 1:midSize
    massFract{j} = strcat(['m' num2str(j-1)]);
end
tBox = t(ones(midSize,1),:);

h1 = plot(tBox',dataMet{toDoF,2}');
for j = 1:midSize
    set(h1(j),'Color',hsvMap(j,:));
end
hold on
for j = 1:midSize
    for k = 1:numel(t)
        lowB = dataMet{toDoF,2}(j,k) - 2*dataMet{toDoF,3}(j,k);
        highB = dataMet{toDoF,2}(j,k) + 2*dataMet{toDoF,3}(j,k);
        plot([t(k) t(k)],[lowB highB], '-','Color',hsvMap(j,:));
    end
end
hold off
title(dataMet{toDoF,1});
ylabel('conc.')
xlabel('time')
legend(massFract(1:midSize));

% 
% subplot(1,2,2);
% h2 = plot(tBox',dataMetB{toDoF,2}');
% for j = 1:midSize
%     set(h2(j),'Color',hsvMap(j,:));
% end
% hold on
% for j = 1:midSize
%     for k = 1:numel(t)
%         lowB = dataMetB{toDoF,2}(j,k) - 2*dataMetB{toDoF,3}(j,k);
%         highB = dataMetB{toDoF,2}(j,k) + 2*dataMetB{toDoF,3}(j,k);
%         plot([t(k) t(k)],[lowB highB], '-','Color',hsvMap(j,:));
%     end
% end
% hold off
% title(strcat([dataMetB{toDoF,1} '\_basal']));
% legend(massFract);
% xlabel('time (min)')
% 
% h1p = get(h1(1),'Parent');
% h2p = get(h2(1),'Parent');
% h1py = get(h1p,'YLim');
% h2py = get(h2p,'YLim');
% if h1py(2) > h2py(2)
%     set(h2p,'YLim',h1py);
% else
%     set(h1p,'YLim',h2py);
% end

drawFigFH = @(sourceX,callbackdataX)drawFig0(sourceX,callbackdataX,dataMet,t);

PBn = uicontrol('Style', 'pushbutton','String','next','Position',[80 20 60 20],...
    'Callback',drawFigFH,'Value',toDoF,'Max',toDoF);
PBb = uicontrol('Style', 'pushbutton','String','back','Position',[20 20 60 20],...
    'Callback',drawFigFH,'Value',toDoF,'Max',toDoF);