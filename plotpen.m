function plotpen(DATA, x, y, varargin)
%
%function plotpen(map, x, y)
%
%plotpen plots the RF locations for all cells recorded from
%one grid co-ordinate, specified by x,y.
%map is an array of structures describing the RFs, as returned by
%loadpenmap


DEFINITIONS;
map = DATA.map;
RF= 1; GRID = 2;
type = RF;
colorbyrf = 0;
j = 1;
colors = mycolors;
while j < nargin -2
    if(strncmpi(varargin{j},'grid',4))
        type = GRID;
    elseif(strncmpi(varargin{j},'colors',4))
        j = j+1;
        colors = varargin{j};
        colorbyrf = 1;
    end
    j = j+1;
end



pens = map.pen;
idx = find(pens(:,2) == x & pens(:,3) == y);

if isempty(idx)
  fprintf('There are no cells at %.1f, %.1f\n',x,y);
  return;
end



if type == RF
    fprintf('Pen at %.1f %.1f\n',x,y);
    rfs = map.rf(idx,:);
    if DATA.plot.popuprf
        plot(rfs(:,1),rfs(:,2),'o');
    end
    hold on;
    for j = 1:size(rfs,1)
        id = idx(j);
        if DATA.selected(idx(j))
                area = map.area(idx(j));
            if colorbyrf
                c = colors{j};
            else
                c = colors{area};
            end
            if DATA.plot.popuprf
                plotrf(rfs(j,:),c);
                h = text(rfs(j,1),rfs(j,2),[map.cellname{idx(j)} ' ' num2str(map.depth(idx(j)))]);
                set(h,'color',c);
            end
            rfsize = sqrt(map.rf(id,3)^2 + map.rf(id,4)^2);
            fprintf('%s %d %d %s %s %.2f %.2f %.0f sz=%.2f\n',map.cellname{idx(j)},...
                map.pen(idx(j,1)),map.area(idx(j)),areanames{map.area(idx(j))},map.datestr{idx(j)},rfs(j,1),rfs(j,2),map.depth(idx(j)),rfsize);

        end
    end
elseif type == GRID
    fprintf('Pen at %.1f %.1f\n',x,y);
    for j = 1:length(idx)
        if DATA.selected(idx(j))
            fprintf('%s %d %s %s %.1f %.1f or %.0f %.1f\n',map.cellname{idx(j)},...
                map.pen(idx(j),1),areanames{map.area(idx(j))},map.datestr{idx(j)},...
                map.rf(idx(j),1),map.rf(idx(j),2),map.rf(idx(j),5),map.depth(idx(j)));
        end
    end
end



function plotrf(rf,col)

x = [-rf(3)  rf(3) rf(3) -rf(3) -rf(3)];
y = [-rf(4)  -rf(4) rf(4) rf(4) -rf(4)];
or = rf(5) * pi/180;
xp = x .* cos(or) + y .* sin(or);
yp = y .* cos(or) - x .* sin(or);
xp = xp+rf(1);
yp = yp+rf(2);

plot(xp,yp,'color',col);
