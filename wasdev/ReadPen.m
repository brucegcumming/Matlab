function pen = ReadPen(file, varargin)

pen = [];
if ~exist(file)
    mycprintf('blue','%s Does not exist\n',file);
    return;
end
fin = fopen(file,'r');
if fin < 0
    fprintf('Can''t Read %s\n',file);
    return;
end

depths = [-2000];
times = [0];
entertime = [];
enterdepth = 0;
DepthOffset = 0;
cellctr = 0;
rwctr = 0;
ncomm = 0;
nrf = 0;
pdate = '??';
noplot = 0;
name = {};
cellid = [];
totaltrials = 0;
j = 1;
while j <= nargin -1
    if strncmpi(varargin{j},'noplot',3)
        noplot = 1;
    end
    j = j+1;
end
toff = 0;
MapDefs;
pen.comments = {};
pen.arealines = [];
pen.missed = 0;
nlines = 0;
nexpts = [];
%toff not tracked properly yet.  Try just using datenum at open + time
%? may need to check for reopen of log without change in t...

while ~feof(fin)
    inlin = fgets(fin);
    nlines = nlines+1;
    if strncmpi(inlin,'ed ',3)
        d = sscanf(inlin,'ed %f at %f');
        if length(d) > 1
        if  length(times) & (d(2)/3600 + toff) < times(end) 
            toff = toff + times(end) - (d(2)/3600);
        end
        times = [times (toff +(d(2)/3600))];
        depths = [depths d(1)+DepthOffset];
        nexpts(length(times)) = 0;
        end
    elseif strfind(inlin,'Miss')
        if strfind(inlin,'Missed MT');
            if pen.missed ~= 2
                pen.missed = MT;
            end
        elseif strfind(inlin,'Missed STS');
            pen.missed = STS;
        else
            pen.missed = VUNKNOWN;
        end
        ncomm = ncomm+1;
        pen.comments{ncomm} = inlin;
        pen.cmtype(ncomm) = 7;
        pen.cmtime(ncomm) = length(times);
    elseif strfind(inlin,'cm=:')
        ncomm = ncomm +1;
        d = strfind(inlin,'cm=:');
        pen.comments{ncomm} = inlin(d+4:end);
        pen.cmtime(ncomm) = length(times);
        pen.cmtimestr(ncomm,:) = inlin(2:10);
        pen.cmtype(ncomm) = 1;
    elseif strfind(inlin,'cm=rf')
        if ~nrf | length(times) > pen.rfs(nrf).time
            nrf = nrf+1;
        end
        [rf] = sscanf(inlin,'cm=rf %f,%f:%fx%f,%fdeg');
        if length(rf) > 4
        pen.rfs(nrf).pos = rf(1:5);
        pen.rfs(nrf).time = length(times);
        else
            pen.rfs(nrf).pos = '';
            pen.rfs(nrf).time = NaN;    
        end
    elseif strfind(inlin,'cm=')
        ncomm = ncomm +1;
        d = strfind(inlin,'cm=');
        pen.comments{ncomm} = inlin(d+3:end);
        pen.cmtime(ncomm) = length(times);
        pen.cmtimestr(ncomm,:) = inlin(2:10);
        pen.cmtype(ncomm) = 2;
            
    elseif strncmpi(inlin,'Opened',6)
        pdate = inlin(8:17);
        pen.datenum = datenum(inlin(12:31));
        p = strfind(inlin,' pen');
        if p
            a = sscanf(inlin(p(1):end),' pen %d %f,%f');
            if length(a) > 2
            pen.num = a(1);
            pen.pos = a(2:3);
            else
                fprintf('%s:Cant parse %s',file,inlin);
            end
        end
    elseif strfind(inlin,'Matter') | strfind(inlin,'WM') | strfind(inlin,'GM')
        ncomm = ncomm+1;
        pen.comments{ncomm} = inlin;
        pen.cmtype(ncomm) = 5;
        pen.cmtime(ncomm) = length(times);
    elseif strfind(inlin,'Boundary')
        strs = split(inlin);
        pen.boundary(1).depth = sscanf(strs{2},'%f'); 
        pen.boundary(1).above = strs{3};
        pen.boundary(1).below = strs{4};
    elseif strfind(inlin,'VisualArea')
        ncomm = ncomm+1;
        pen.comments{ncomm} = inlin;
        pen.cmtype(ncomm) = 6;
        pen.cmtime(ncomm) = length(times);
        pen.arealines(end+1) = ncomm;
    elseif strfind(inlin,'Impedance')
        ncomm = ncomm+1;
        pen.comments{ncomm} = inlin;
        pen.cmtype(ncomm) = 4;
        pen.cmtime(ncomm) = length(times);
    elseif strfind(inlin,'File')
        cellctr = cellctr+1;
        name{cellctr} = splitpath(inlin);
        cellid(cellctr) = length(times);
        
    elseif strncmpi(inlin,'StartDepth ',8)
        [d, entertime] = sscanf(inlin,'StartDepth %f at %f');
        if isempty(d)
            d = sscanf(inlin,'StartDepth %f');
        else
            id = findstr(inlin,'at ');
            entertime = inlin(id+4:end-1);
            enterdepth = d(1);
        end
        if isempty(d)
            depths = [depths 0];
            times =[times 0];
        else
            depths = [depths d(1)];
            times = [times times(end)];
        end
        nexpts(length(times)) = 0;
    elseif strncmpi(inlin,'Expt',4)
        if length(times) > 1
            nexpts(length(times)) = nexpts(length(times))+1;
        end
    elseif strncmpi(inlin,'Rewards ',7)
        d = sscanf(inlin,'Rewards %f of %f');
        rwctr = rwctr+1;
        goodtrials(rwctr) = d(1);
        if length(d) > 1
        totaltrials(rwctr) = d(2);
        end
        if length(times)
            rwid(rwctr) = length(times);
        else
            rwid(rwctr) = 1;
        end
     
    elseif strfind(inlin,'File')
    end
    if isfield(pen,'datenum')
        pen.date(length(times)) = pen.datenum+(times(end)-toff)./(24);
    end
end

%since times are taken from depth changes etc, not every line
%time of file should match time of comment
if ~isempty(pen.arealines)
    al = pen.cmtime(pen.arealines);
    for j = 1:length(cellid)
        pen.visualarea{j} = 'Vd';
        id = find(al <= cellid(j));
        if ~isempty(id)
            s = pen.comments{pen.arealines(id(end))};
            xid = strfind(s,'VisualArea');
            if ~isempty(xid);
                pen.visualarea{j} = regexprep(s(xid+11:end),' .*','');
            end
        end
    end
end
pen.times = times;
pen.nexpts = nexpts;
pen.depths = depths;
pen.entertime = entertime;
pen.enterdepth = enterdepth;
pen.totaltrials = totaltrials;
pen.files = name;
pen.filetimes = cellid;
fclose(fin);
if noplot
    return;
end

hold off;
plot(times(2:end),depths(2:end));
for j = 1:cellctr
    text(times(cellid(j)),depths(cellid(j)),name{j},'Rotation',90);
end
hold on;
plot(times(rwid),totaltrials,'r');
scale = max(depths)/(1+max(nexpts)*2);
lax = gca;
legend('depth','#trials');
rax = AddRplot(lax,times(2:end),cumsum(nexpts(2:end)),'g');
legend('#expts');
title(sprintf('%s %s %s',file,pdate,entertime));
pen.axes = [lax rax];