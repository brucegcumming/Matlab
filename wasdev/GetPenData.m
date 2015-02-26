function [area, depth, pen] = GetPenData(name, varargin)

global bgcfileprefix;
global cellmap;
reload = 0;
area = [];
depth = [];
pen = [];

j = 1;
while j < nargin
   if strncmpi('reload',varargin{j},4)
       reload = 1;
   end
   j = j+1;
end

if isempty(cellmap) | reload
    dfile = '/bgc/bgc/anal/dufus/dufv1.fixtab';
   if exist(dfile,'file') 
    [drx, dry, drw, drh, dro, dnames, dids, dpx, dpy, darea, ddepth, ddates] = textread(dfile,'%f%f%f%f%f%s%d%f%f%s%d%s%*[^\n]');
   else
       dnames = [];
   end
    rfile = '/bgc/bgc/anal/rufus/rufv1.fixtab';
   if exist(rfile,'file') 
    [rrx, rry, rrw, rrh, rro, rnames, rids, rpx, rpy, rarea, rdepth, rdates] = textread(rfile,'%f%f%f%f%f%s%d%f%f%s%d%s%*[^\n]');
   else
       rnames = [];
   end
   if length(rnames) & length(dnames)
    cellmap.names = [dnames; rnames];
    cellmap.area = [darea; rarea];
    cellmap.depth = [ddepth; rdepth];
    cellmap.pens = [dids; rids];
    cellmap.penx = [dpx; rpx];
    cellmap.peny = [dpy; rpy];
    for j = 1:length(cellmap.area)
        aid = strmatch(cellmap.area(j),{'V1', 'V2','Vd','??'});
        if isempty(aid) 
            cellmap.areaid(j) = 4;
        elseif  aid == 3
            cellmap.areaid(j) = 1;
        else
            cellmap.areaid(j) = aid;
        end
    end
   end
end


name = splitpath(name,'cellpref');
if isempty(cellmap)
    id = [];
else
    id = strmatch(name, cellmap.names);
end
if isempty(id)
    area = NaN;
    depth = NaN;
else
    if(length(id) > 1)
        fprintf('%s %d entries:',name,length(id));
        fprintf('%s\n',cellmap.names{id});
    end
    area = cellmap.areaid(id(1));
    depth = cellmap.depth(id(1));
    pen.pen = cellmap.pens(id(1));
    pen.x = cellmap.penx(id(1));
    pen.y = cellmap.peny(id(1));
end