function [fits, details]  = CalcPostSaccGain(Expt,varargin)

figtag = 'Postsaccade';

nmin = 10;
sumdur = 500;
smw = 0;
nameprefix = [];
playsafe = 0;
tvals = -1000:100:1500; 
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'dur',3)
        j = j+1;
        sumdur = varargin{j};
    elseif strncmpi(varargin{j},'nmin',4)
        j = j+1;
        nmin = varargin{j};
    elseif strncmpi(varargin{j},'prefix',4)
        j = j+1;
        nameprefix = varargin{j};
    elseif strncmpi(varargin{j},'smooth',4)
        j = j+1;
        smw = varargin{j};
    elseif strncmpi(varargin{j},'tvals',4)
        j = j+1;
        tvals = varargin{j};
    end
    j = j+1;
end

if iscell(Expt)
    if isfield(Expt{1},'varratio')
        PlotFits(Expt,varargin{:});
        return;
    end
    for j = 1:length(Expt)
    if ischar(Expt{j})
        E = LoadExpt(Expt{j});
    else
        E = Expt{j};
    end
    [fits{j}, details{j}] = CheckPostSaccGain(E, varargin{:});
    end
    return;
end

if ischar(Expt) && isdir(Expt)
    nameprefix = Expt;
    a = dir([Expt '/*OPRC.cell*']);
    b = dir([Expt '/*OPRC.mu*']);
    Expt = cat(1,a,b);
end

if isstruct(Expt) && isfield(Expt,'isdir')  %a directory listing
    for j = 1:length(Expt)
        fprintf('Processing %s\n',Expt(j).name);
        E = LoadExpt([nameprefix '/' Expt(j).name]);
        [fits{j}, details{j}] = CalcPostSaccGain(E, varargin{:});
        details{j}.fits = fits{j};
        if playsafe
        save('PostSaccFits.mat','fits','details');
        end
    end
    fits = details;
    return; 
end
SpkDefs;

if smw > 0
    smargs = {'box' smw};
else
    smargs = {'fbox'};
end
details.name = Expt.Header.loadname;
details.cellnumber = Expt.Header.cellnumber;
details.probe = Expt.Header.probe;
[rc, a] = PlotRevCorAny(Expt,'collapse','ce','makeall',smargs{:});
details.varratio = rc.varratio;
[a,b,c] = mksacsdf(Expt.Trials,106,rc.times,'MinSize',0.5,smargs{:});
details.sacsdf = a;
[a,b,c] = PlotRC(rc,'var');
details.meanrate = mean(rc.y(:));
details.var = a;
id = find(~isnan(rc.sdfs.n));
[tdiff, tid] = min(abs(b-rc.delaysamples));
ib = find(rc.sdfs.extraval == IBLANK);
details.resp = [rc.y(id,1,tid)' rc.sdfs.extras{ib}.sdf(b)];
bestt = rc.times(b);
details.xvals = cat(1,rc.x(:,1), IBLANK);
for d = 1:length(sumdur)
    
for j = 1:length(tvals)
    t = tvals(j); 
    Expt = FillTrials(Expt,'saccades',[1 t sumdur(d)]);
    rca = PlotRevCorAny(Expt,'collapse','ce','slices',bestt,'exp2','PostSacc','yval',0,smargs{:},'nmin',nmin,'makeall');
    rcb = PlotRevCorAny(Expt,'collapse','ce','slices',bestt,'exp2','PostSacc','yval',1,smargs{:},'nmin',nmin,'makeall');
    [mresp, mid] = max(rca.y);
    ib = find(rca.sdfs.extraval == IBLANK);
    asdf(1,:) = rca.sdfs.extras{ib}.sdf;
    bsdf(1,:) = rcb.sdfs.extras{ib}.sdf;
    asdf(2,:) = rca.sdfs.s{mid};
    bsdf(2,:) = rcb.sdfs.s{mid};
    id = find(~isnan(rca.sdfs.n+rcb.sdfs.n));
    x = [rca.y(id)' rca.extray(1)];
    y = [rcb.y(id)' rcb.extray(1)];
    fits(j,:) = polyfit(x,y,1);
    details.sacresps(j,:) = y;
    details.n(j,:) = [sum(rca.sdfs.n) sum(rcb.sdfs.n)];
end
details.fits = fits; 
details.ctime = now;
details.name = Expt.Header.loadname;
details.nmin = nmin;
details.smw = smw;
details.sumdur = sumdur./10;
details.delay = bestt/10;


gains(d,:) = fits(:,1);
offsets(d,:) = fits(:,2);
details.tvals= tvals./10;
GetFigure(figtag);
subplot(1,1,1);
hold off;
if length(tvals) > 1
[ax, ha, hb] = plotyy(details.tvals,fits(:,1),details.tvals,fits(:,2));
set(ha,'marker','o');
set(hb,'marker','o');
else
    hold off; 
    plot(rca.times./10,asdf(1,:)); hold on; plot(rcb.times./10,bsdf(1,:),'r');
    plot(rca.times./10,asdf(2,:),'--');  plot(rcb.times./10,bsdf(2,:),'r--');
    plot([bestt bestt]/10,get(gca,'ylim'),'k');
end
end
if length(sumdur) > 1
    subplot(2,1,1);
    imagesc(gains);
    colorbar;
    subplot(2,1,1);
    imagesc(offsets);
    colorbar;
end

function PlotFits(fits, varargin)


DATA.minvar = 10;
DATA.minrate = 1;
plottype = 'twoimage';
j = 1;
while j < length(varargin)
    if strncmpi(varargin{j},'minvar',4)
        j = j+1;
        DATA.minvar = varargin{j};
    elseif strncmpi(varargin{j},'minrate',4)
        j = j+1;
        DATA.minrate = varargin{j};
    elseif strncmpi(varargin{j},'plot',4)
        j = j+1;
        plottype = varargin{j};
    end
    j = j+1;
end
for j = 1:length(fits)
    if ~isfield(fits{j},'delay')
    fits{j}.delay = 60;
    end
    depths(j) = fits{j}.probe;
end

[a,did] = sort(depths);
DATA.depthid = did;

DATA.tag.mainplot = 'PostSaccade';
DATA.tag.onefit = 'OneSaccFit';
GetFigure(DATA.tag.mainplot);
if strcmp(plottype,'varim')
    hold off;
    n = 0;
    for j = 1:length(fits)
        k = did(j);
        v = fits{k}.var;
        if fits{k}.varratio(1) > DATA.minvar && fits{k}.meanrate > DATA.minrate
            n = n+1;
            DATA.vars(n,:) = v./mean(v);
            DATA.imrow(n) = k;
        end
    end
    imagesc(DATA.vars,'ButtonDownFcn',{@HitImage, 'var'});
elseif strcmp(plottype,'var')
    clf;
    hold off;
    for j = 1:length(fits)
        h = plot(fits{j}.varratio(1),fits{j}.varratio(2),'o','ButtonDownFcn',{@HitPoint, j});
        if strfind(fits{j}.name,'.cell')
            set(h,'color','r');
        end
        hold on;
    end
elseif strcmp(plottype,'varmean')
    clf;
    hold off;
    for j = 1:length(fits)
        plot(fits{j}.varratio(1),fits{j}.meanrate,'o','ButtonDownFcn',{@HitPoint, j});
        hold on;
    end
else   
n = 0;
for j = 1:length(fits)
    k = did(j);
    if fits{k}.varratio(1) > DATA.minvar && fits{k}.meanrate > DATA.minrate 
        n = n+1;
        g = fits{did(j)}.fits(:,1);
        o = fits{did(j)}.fits(:,2);
        if isnan(sum(g))
            id = find(isnan(g));
            gid = find(~isnan(g));
            g(id) = interp1(gid,g(gid),id);
            o(id) = interp1(gid,o(gid),id);
        end
        DATA.fitgains(n,:) = g./mean(g);
        DATA.fitoffsets(n,:) = o./fits{k}.meanrate;
        DATA.imrow(n) = k;
    end
end


subplot(2,1,1);
h = imagesc(DATA.fitgains);
set(h,'buttondownfcn',{@HitImage, 'gains'});
colorbar;
subplot(2,1,2);
h = imagesc(DATA.fitoffsets);
colorbar;
set(h,'buttondownfcn',{@HitImage, 'offsets'});
GetFigure(DATA.tag.onefit);
clf;
plot(fits{1}.tvals, mean(DATA.fitgains));
GetFigure(DATA.tag.mainplot);
end
DATA.fits = fits;
DATA.toplevel = gcf;
set(DATA.toplevel,'UserData',DATA);


function HitPoint(a,b,cid)

DATA = GetDataFromFig(a);
fit = DATA.fits{cid};
Expt = Loadexpt(fit.name);
[a,b] = PlotRevCorAny(Expt,'collapse','ce','slices',fit.delay*10,'nmin',fit.nmin,'makeall','fbox');
figure(a.figa);
set(a.legend,'position',[0.8 0.5 0.2 0.5]);
title(sprintf('%s P%.1f',CellName(fit.name),fit.probe));

function s = CellName(name)

s = strrep(name,'.mat','');
s = regexprep(s,'.*\.','');



function HitImage(a,b,type)

DATA = GetDataFromFig(a);
SpkDefs;

pos = get(gca,'currentpoint');
j = round(pos(1,2));
k = round(pos(1,1));
shiftalt = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
fit = DATA.fits{DATA.imrow((j))};
GetFigure(DATA.tag.onefit);
subplot(2,2,1);
plotyy(fit.tvals,fit.fits(:,1),fit.tvals,fit.fits(:,2));
[a,b] = fileparts(fit.name);
title(sprintf('%s P%.1f',b,fit.probe));
subplot(2,2,2);
hold off;
id = find(fit.xvals > -1000);
bid = find(fit.xvals == IBLANK);
plot(fit.xvals(id),mean(fit.sacresps(:,id)));
hold on;
plot(fit.xvals(id),fit.sacresps(k,id),'m');
plot(fit.xvals(id),fit.resp(id),'r');
if length(bid) ==1
    plot(minmax(fit.xvals(id)),[fit.resp(bid) fit.resp(bid)],'r--');
    plot(minmax(fit.xvals(id)),[fit.sacresps(k,bid) fit.sacresps(k,bid)],'m--');
end    
subplot(2,2,3);
imagesc(fit.sacresps);
subplot(2,2,4);
plot(fit.resp, fit.sacresps(k,:),'o');
refline(1);
Expt = Loadexpt(fit.name);
if strcmp(type,'var')
    [a,b] = PlotRevCorAny(Expt,'collapse','ce','slices',fit.delay*10,'fbox','nmin',fit.nmin,'makeall');
else
    Expt = FillTrials(Expt,'saccades',[1 fit.tvals(k)*10 fit.sumdur*10]);
    [a,b] = PlotRevCorAny(Expt,'collapse','ce','slices',fit.delay*10,'exp2','PostSacc','fbox','nmin',fit.nmin,'makeall');
end
figure(a.figa);
set(a.legend,'position',[0.8 0.5 0.2 0.5]);
plot(a.times./10,fit.sacsdf,'k-','linewidth',2);
title(sprintf('%s P%.1f Var %.1f,%.1f',CellName(fit.name),fit.probe,fit.varratio(1),fit.varratio(2)));

