function PlotAllExpts(a,b, varargin)
if isfield(a, 'state') && isfield(a,'AllData')
DATA = a;
else
DATA = GetDataFromFig(a);
end
figlabela = 'AllExpts';
DATA.firsttrial = 0;
DATA.lasttrial = 0;
j = 1;
while j <= length(varargin)
if strncmpi(varargin{j},'figlabel',4)
j = j+1;
figlabela = varargin{j};
end
j = j+1;
end
nex = length(DATA.Expts);
[nr, nc] = Nsubplots(nex);
pw = 0.9./nc;
ph = (1-0.02.*nr)./nr;
hstep = 0.02;
vstep = 0.02;
set(DATA.toplevel,'Name','Busy......');
drawnow;

DATA.ptsize = 1;
[f, isnew]  = GetFigure(figlabela);
if isnew  
hm = uimenu(f,'Label','Spool','Tag','XClusterMenu');
uimenu(hm,'Label','1/10','Callback',{@cmb.SpoolAllExpts, DATA.toplevel, 1});
uimenu(hm,'Label','1/20','Callback',{@cmb.SpoolAllExpts, DATA.toplevel, 2});
uimenu(hm,'Label','1/50','Callback',{@cmb.SpoolAllExpts, DATA.toplevel, 3});
end


if strncmp(DATA.filetype,'Grid',4)
return;
end

tstart = now;
for xi = 1:nex
y = vstep + floor((xi-1)/nc) .* (ph+vstep);
x = hstep + (mod(xi-1, nc))./(nc * 1.02);
if isfield(DATA.Expts{xi}.gui,'spks')  %force, in case reloaded clusters
[DATA, ispk] = SetExptSpikes(DATA,xi, 0,'useexpall');
elseif 1
[DATA, ispk] = SetExptSpikes(DATA,xi, 0,'useexp');
end
if ~isempty(ispk)
DATA.nclusters = cmb.CountClusters(DATA.cluster);
set(0,'CurrentFigure',f);
subplot('Position' ,[x y pw ph]);
hold off;
cmb.DrawXYPlot(DATA, DATA.Expts{xi}.gui.spks,'usefigure');
xl = get(gca,'Xlim');
yl = get(gca,'ylim');
if cmb.iscluster(DATA.cluster,1,DATA.probe) & isfield(DATA.cluster{1,DATA.probe},'dprime')
dp = DATA.cluster{1,DATA.probe}.dprime;
text(xl(1)+diff(xl)/5,yl(2)-diff(yl)/10,sprintf('%d:%.1f',xi,dp),'color','r');
else
text(xl(1)+diff(xl)/5,yl(2)-diff(yl)/10,sprintf('%d ed=%.1f',xi,GetEval(DATA.Expts{xi},'ed')),'color','r');
end
%   drawnow;
xlabel('');
ylabel('');
title('');
end
end
mytoc(tstart);
cmb.NotBusy(DATA);
set(DATA.toplevel,'UserData',DATA);

