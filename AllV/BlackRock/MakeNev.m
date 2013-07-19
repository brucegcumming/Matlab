function out = MakeNev(name, varargin)
%GUI for building Nev files from.ns5 files

TOPTAG = 'Ns2Nev';
it = findobj('Tag',TOPTAG);
prctile = 0;
sdth = 0;
setstate = [];
j = 1;
while j <= length(varargin)
    if ischar(varargin{j})
    if strncmpi(varargin{j},'prc',3)
        j = j+1;
        prctile = varargin{j};
    elseif strncmpi(varargin{j},'sdth',4)
        j = j+1;
        sdth = varargin{j};
    elseif strncmpi(varargin{j},'writenev',8)
        setstate.writenev = 1;
    end
    end
    j = j+1;
end
if isempty(it)  %new GUI/

DATA.toplevel= GetFigure(TOPTAG);
DATA.tags.top = TOPTAG;
DATA.tags.spikev = [TOPTAG 'SpikeV'];
DATA.nprobes = 1;
DATA.mouse.down = 0;
DATA.lowpass = Gauss(25,-200:200);
f = fields(setstate);
for j = 1:length(f)
    DATA.state.(f{j}) = setstate.(f{j});
end
if isstruct(name) && isfield(name,'MetaTags')
    nsx = name;
elseif isdir(name)
    DATA = BuildGUI(DATA);
    d = dir(name);
    nf = 0;
    for j = 1:length(d)
        if strfind(d(j).name,'.ns5')
            nf = nf+1;
            filenames{nf} = [name '/' d(j).name];
            out{nf} = MakeNev(filenames{nf},varargin{:});
            out{nf}.name = filenames{nf};
        end
    end
    tsum = 0;
    esum = 0;
    for j = 1:nf
        tsum = tsum +out{nf}.took;
        esum = esum +out{nf}.duration;
    end
    fprintf('Total time for directory %.2f, took %.2f\n',esum,tsum)
    return;
else
    DATA = BuildGUI(DATA);
    [DATA, out] = ProcessFile(DATA, name, varargin{:});
    return;
end

DATA.nsx = nsx;
for j = 1:size(nsx.Data,1)
   v = conv(nsx.Data(j,:),DATA.lowpass);
   v = v(201:end-200);
   DATA.nsx.Data(j,:) = DATA.nsx.Data(j,:)-v; 
end
DATA.nprobes = size(nsx.Data,1);
DATA.npts = size(nsx.Data,2);
DATA = BuildGUI(DATA);
set(DATA.toplevel,'UserData',DATA);
set(gcf, 'WindowButtonDownFcn',@ButtonPressed);
set(gcf, 'WindowButtonUpFcn',@ButtonReleased);
set(gcf, 'WindowButtonMotionFcn',@ButtonDragged);
else
    DATA = get(it,'UserData');
    if strncmpi(name,'peakf',5)
        DATA.peakf = varargin{1}; %actually period
        if length(varargin) > 1
            DATA.filtersd = varargin{2}*DATA.peakf;
        end
        DATA = TestPlot(DATA,'prct',1);
        set(DATA.toplevel,'UserData',DATA);
    elseif strncmpi(name,'prc',3)
        prctile = varargin{1};
    elseif strncmpi(name,'plotpcs',7)
        PlotSpikes(DATA,'plotpcs',varargin{:});
    elseif exist(name,'file')
        [DATA,out] = ProcessFile(DATA, name, varargin{:});
   end
    if prctile > 0
        ApplyAndPlot(DATA,'prct',prctile);
    end
end


function [DATA, out] = ProcessFile(DATA, name, varargin)

hdr = openNSx(name);
DATA.filename = name;
DATA.nprobes = 1;
DATA.nev.MetaTags = hdr.MetaTags;
tstart = now;
for j = 1:length(hdr.MetaTags.ChannelID)
    e = sprintf('e:%d',j);
    ts = now;
    fprintf('Reading Channel %d...',j);
    DATA = ReadProbeFromNs5(DATA,j,varargin{:});
    fprintf('took %.3f %d Spikes, %.1fHz\n',mytoc(ts),size(DATA.spks,1),size(DATA.spks,1)./DATA.duration);
    if DATA.state.writenev
    WriteNev(DATA);
    end
    out.probe(j).nspk = size(DATA.spks,1);
    out.probe(j).sd = std(DATA.nsx.Data(DATA.probe,:));
end
out.duration = DATA.duration;
out.took = mytoc(tstart);
fprintf('%s %d channels duration %.2f took %.1f\n',DATA.filename,j,DATA.duration,mytoc(tstart));
return;


function DATA = ReadProbeFromNs5(dat, probe, varargin)

    e = sprintf('e:%d',probe);
    dat.nsx = openNSx(dat.filename,'read',e);
    v = dat.nsx.Data(1,:);
    npts = length(v);
    ml = floor(ceil(npts/500)).*500;
    v(length(v)+1:ml) = 0;
    mains = mean(reshape(v, 500,ml/500),2)';
    v = v - repmat(mains,1,ml/500);
    v = v(1:npts);
    lf = conv(v,dat.lowpass);
    lf = lf(201:end-200);
    dat.nsx.Data(1,:) = v-lf;
    dat.nsx.Data(1,1:200) = 0;
    dat.nsx.Data(1,end-200:end) = 0;
    dat.probeid = probe;
    dat.npts = size(dat.nsx.Data,2);
    dat.probe = 1;
    dat.duration = dat.npts./30000;

    DATA = ApplyAndPlot(dat, varargin{:});
    

function WriteNev(DATA, varargin)

oname = strrep(DATA.filename,'.ns5',sprintf('p%d.mat',DATA.probeid));
[a,b,c] = fileparts(DATA.filename);
oname = [a '/GridSpikes/' b sprintf('p%d.mat',DATA.probeid)];
nev = DATA.nev;
nev.Data.Spikes.Waveform = DATA.spks;
nev.Data.Spikes.TimeStamp = DATA.TimeStamps';
nev.Data.Spikes.codes = zeros(size(DATA.spks,1),1);
nev.probeid = DATA.probeid;
save(oname,'nev','-v7.3');


function DATA = ApplyAndPlot(DATA, varargin)

plottype = 2;
j = 1;
type = 0;
while j <= length(varargin)
    if strncmpi(varargin{j},'pca',3)
        plottype = 2;
    elseif strncmpi(varargin{j},'prct',4)
        j = j+1;
        DATA.lowth = prctile(DATA.nsx.Data(DATA.probe,:),varargin{j})
    elseif strncmpi(varargin{j},'sdth',4)
        j = j+1;
        DATA.lowth = varargin{j};
        type = 1;
    end
    j = j+1;
end
[DATA.spks, DATA.TimeStamps] = ApplyThreshold(DATA.nsx.Data(DATA.probe,:),DATA.lowth, 90, 30, type);

PlotSpikes(DATA);

function PlotSpikes(DATA, varargin)

plottype = 2;
pcplot = [1 2];
plotwave = 1;
plotdensity = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'density',3)
        plotdensity = 1;
    elseif strncmpi(varargin{j},'pca',3)
        plottype = 2;
    elseif strncmpi(varargin{j},'plotpcs',7)
        plottype = 2;
        j = j+1;
        pcplot = varargin{j};
        plotwave = 0;
    elseif strncmpi(varargin{j},'prct',4)
        j = j+1;
        DATA.lowth = prctile(DATA.nsx.Data(DATA.probe,:),varargin{j})
    end
    j = j+1;
end

if plotwave
GetFigure(DATA.tags.spikev,'front');
plot(DATA.spks','color',[0.5 0.5 0.5]);
title(sprintf('%d Spikes %.2fHz',size(DATA.spks,1),size(DATA.spks,1)./DATA.duration));
end
GetFigure('XYPlot','front');
ptsz = 1;
if size(DATA.spks,1) < 1000
    ptsz = 2;
end
if plottype == 3
    hist(v)
elseif plottype == 2
[evec,eVal] = eig(cov(DATA.spks));
pcs = DATA.spks * evec;
plot(pcs(:,end-pcplot(1)+1),pcs(:,end-pcplot(2)+1),'.','markersize',ptsz);
else
dVdt=diff(DATA.spks,[],2);
E = sum(dVdt.^2,2);
V = var(DATA.spks,[],2);
plot(E,V./E,'.','markersize',1);
end

if plotdensity
DensityPlot(pcs(:,end-pcplot(1)+1),pcs(:,end-pcplot(2)+1),'sd',[2 2],'nbins',[100 100]);
end
set(DATA.toplevel,'UserData',DATA);


function DATA = TestPlot(DATA, varargin)

if DATA.probe == 8
   w = 200;
   x = -w:w;
   DATA.lowpass = gauss(DATA.filtersd,-w:w) .* cos(x * 2 * pi *1/DATA.peakf);
%   DATA.lowpass = -diff(diff(gauss(DATA.filtersd,-(w+1):(w+1))));
   
   v = conv(DATA.nsx.Data(7,:),DATA.lowpass);
   DATA.nsx.Data(8,:) = v(w+1:end-w); 
   GetFigure(DATA.toplevel);
   hold off;
   plot(DATA.nsx.Data(DATA.probe,:));
   DATA = ApplyAndPlot(DATA,varargin{:});
   GetFigure(DATA.tag.spikev,'front');
   hold on;
   plot(DATA.lowpass(w:w+89).*max(DATA.spks(:))/max(DATA.lowpass));
   hold off;
end

function DATA = BuildGUI(DATA)

bp(1) = 0.02;
bp(2) = 0.01;
bp(3) = 0.1;
bp(4) = 0.05;
uicontrol(DATA.toplevel,'style','pop','string',num2str([1:DATA.nprobes]'),'Units','Norm','Position',bp,'Tag','Probeid',...
        'Callback',@NewProbe);
    set(DATA.toplevel,'UserData',DATA);

    
    
function NewProbe(a,b)

DATA = GetDataFromFig(a);
DATA.probe = get(a,'value');
DATA.duration = DATA.npts./30000;
GetFigure(DATA.toplevel);
hold off;
plot(DATA.nsx.Data(DATA.probe,:));
set(DATA.toplevel,'UserData',DATA);

function ButtonPressed(src, data)

DATA = GetDataFromFig(src);
start = get(gca,'CurrentPoint');
DATA.npts = size(DATA.nsx.Data,2);
DATA.tag.spikev = 'SpikeV';
if isfield(DATA,'lowh') && ishandle(DATA.lowh)
    set(DATA.lowh,'Xdata',[1 DATA.npts],'YData',[start(1,2) start(1,2)]);
else
    hold on;
    DATA.lowh = plot([1 DATA.npts],[start(1,2) start(1,2)],'k-');
end
DATA.mouse.down = 1;
set(DATA.toplevel,'UserData',DATA);

function ButtonDragged(src, data)
DATA = GetDataFromFig(src);
if DATA.mouse.down
start = get(gca,'CurrentPoint');
if isfield(DATA,'lowh') && ishandle(DATA.lowh)
    set(DATA.lowh,'Xdata',[1 DATA.npts],'YData',[start(1,2) start(1,2)]);
end
end

function ButtonReleased(src, data)

DATA = GetDataFromFig(src);
start = get(gca,'CurrentPoint');
if isfield(DATA,'lowh') && ishandle(DATA.lowh)
    set(DATA.lowh,'Xdata',[1 DATA.npts],'YData',[start(1,2) start(1,2)]);
    drawnow;
end
DATA.mouse.down = 0;
DATA.lowth = start(1,2);
set(DATA.toplevel,'UserData',DATA);
ApplyAndPlot(DATA);


function [spks, t] = ApplyThreshold(v,th, nsmp, pre, type)

dvdt = diff(v);
post = nsmp-pre;
sgn = diff(sign(dvdt));
if type == 1
    sth = th;
    sd = std(v);
    if length(th) ==1 && th < 0
        th(1) = sd .* sth;
    elseif length(th) == 2
        th(1) = sd .* -sth;
        th(2) = sd .*sth;
    end
end
if th(1) < 0
    t = find(sgn(pre:end-post) > 0 & v(1+pre:end-1-post) < th(1))+pre;

if length(th) > 1 & th(2) > 0
    tu = find(sgn(pre:end-post) < 0 & v(1+pre:end-1-post) > th(2))+pre;
    t = union(tu, t);
    prc = 1./300; %aim for 1 spk/sec
    while size(t,2) == 0
        th = prctile(v,[prc 100-prc]);
        t = find(sgn(pre:end-post) > 0 & v(1+pre:end-1-post) < th(1))+pre;
        tu = find(sgn(pre:end-post) < 0 & v(1+pre:end-1-post) > th(2))+pre;
        t = union(tu, t);
        prc = prc + 1./300;
    end
else
    prc = 1./300; %aim for 1 spk/sec
    while size(t,2) == 0
        th = prctile(v,prc);
        t = find(sgn(pre:end-post) > 0 & v(1+pre:end-1-post) < th(1))+pre;
        prc = prc + 1./300;
    end
end
    
else
t = find(sgn(pre:end-post) < 0 & v(1+pre:end-1-post) > th(1))+pre;
end

%if th < 0
%    t = find(sgn > 0 & v(nsmp:end-nsmp) <= th & v(nsmp-1:end-nsmp-1) > th)' + nsmp-1;
%else
%    t = find(sgn < 0 & v(nsmp:end-nsmp) >= th & v(nsmp-1:end-nsmp-1) < th)' + nsmp-1;
%end
    ts = repmat(t',1,nsmp) + repmat([1:nsmp]-pre,length(t),1);
    spks = v(ts);
    
    