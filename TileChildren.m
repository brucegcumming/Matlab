function TileChildren(parent, varargin)
%resize child windows so taht there are all in the display

closeall = 0;
j = 1;
args = {};
while j <= length(varargin)
    if strncmpi(varargin{j},'closeall',8)
        closeall = 1;
    end
    j = j+1;
end
%find hidden figure too
%figs = findobj('type','figure');
c = allchild(0);
types = get(c,'type');
fid  = find(strcmp('figure',types));
figs = c(fid);
%figs = setdiff(figs,parent);
screens = get(0,'MonitorPositions');
screenpos = sum(screens,1); %w,h are total across screens
if isempty(figs)
    return;
end

for j = 1:length(figs)
    ischild = 0;
    if figs(j) == parent
        figpos(j,:) = get(figs(j),'position');
        ischild  = 2;
    elseif getappdata(figs(j),'ParentFigure') == parent
        figpos(j,:) = get(figs(j),'position');
        ischild = 1;
    end
    if ischild
    newfigpos(j,:) = figpos(j,:);
    if sum(figpos(j,[1 3])) > screenpos(3)
        newfigpos(j,1) = screenpos(3) - figpos(3);
    end
    if sum(figpos(j,[2 4])) > screenpos(4)
        newfigpos(j,2) = screenpos(4) - figpos(3);
    end
    end
end

npix = sum(figpos(:,3).*figpos(:,4));
pixratio = npix./(screenpos(3).*screenpos(4));
for j = 1:length(figs)
    if sum(newfigpos(j,:))
        set(figs(j),'position',newfigpos(j,:));
    end
end
end