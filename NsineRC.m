function [res, ns] = NsineRC(Expt, All, varargin)

%
% NsineRC(Expt, All, delays, ids, ....)
% builds RC plots for data with Nsine stimuli
%

if nargin > 2 & isnumeric(varargin{1}) % delays
    delays = varargin{1};
    j = 2;
    if nargin > 3 & iscell(varargin{2})
        ids = varargin{2};
        j = 3;
    end
else
    delays = [];
    ids = {};
    j = 1;
end
sfid = 0;
ns = [];
plottype = 0;
sfvals =[];
seedoffset =0;
eyes = 0;
MONOC = 5;
DPSF = 6;
logscale = [0 0];

verbose = 0;
xv = [];
if isempty(ids)
    ids = {}; %%needs to be a cell
end

while j < nargin-1
  if strncmpi(varargin{j},'delays',3)
    j = j+1;
    delays = varargin{j};
  elseif strncmpi(varargin{j},'dpsf',3)
      plottype = DPSF;
  elseif strncmpi(varargin{j},'fta',3)
      plottype =2;
  elseif strncmpi(varargin{j},'logx',4)
      logscale(1) = 1;
  elseif strncmpi(varargin{j},'offset',2)
    j = j+1;
    seedoffset = varargin{j};
  elseif strncmpi(varargin{j},'sta',3)
      plottype =1;
  elseif strncmpi(varargin{j},'monoc',3)
    plottype = MONOC;
  elseif strncmpi(varargin{j},'rc',2)
    j = j+1;
    rcres = varargin{j};
    All.rcid = rcres.rcid;
    All.delays = rcres.delays;
  elseif strncmpi(varargin{j},'xv',2)
    j = j+1;
    xv = varargin{j};
  elseif strncmpi(varargin{j},'sf',2)
    j = j+1;
    sfvals = varargin{j};
  elseif strncmpi(varargin{j},'phases',3)
      plottype = 3;
      j = j+1;
      sfid = varargin{j};
  end
  j = j+1;
end

if isempty(delays)
    if isfield(All,'delays') & isfield(All,'rcid')
        ids = All.rcid;
        delays = All.delays;
    else
        delays = [200 250 300 350 400 450 500 550 600 650 700 750 800];
    end
end

if isfield(All,'seeds')
    ns = All;
else
    if isinf(All(3,1))
        ns.seeds = All(1:4:end,1);
        ns.ids = All(2:4:end,1);
        ns.lsfs = All(1:4:end,2:end);
        ns.phases = All(2:4:end,2:end) * pi/180;
        ns.rsfs = All(3:4:end,2:end);
        ns.sfs = ns.lsfs & ns.rsfs;
        ns.dp = All(4:4:end,2:end) * pi/180;
        ns.dpvals = unique(ns.dp(:));
    elseif isnan(All(3,1))
        ns.seeds = All(1:3:end,1);
        ns.ids = All(2:3:end,1);
        ns.sfs = All(1:3:end,2:end);
        ns.phases = All(2:3:end,2:end) * pi/180;
        ns.dp = All(3:3:end,2:end) * pi/180;
        ns.dpvals = unique(ns.dp(:));
    else
        ns.seeds = All(1:2:end,1);
        ns.ids = All(2:2:end,1);
        ns.sfs = All(1:2:end,2:end);
        ns.phases = All(2:2:end,2:end) * pi/180;
    end
    ns.name = splitpath(Expt.Header.Name);
end

if ismember(plottype,[1 2])
%
% N.B dp should be in radians, not degrees
    if isfield(ns,'dp')
        lp = ns.phases - ns.dp/2;
        rp = ns.phases + ns.dp/2;
        coss = ns.sfs .* cos(lp);
        sins = sin(lp) .* ns.sfs;
        if isfield(ns,'rsfs')
            rcoss = ns.rsfs .* cos(rp);
            rsins = sin(rp) .* ns.rsfs;
        else
            rcoss = ns.sfs .* cos(rp);
            rsins = sin(rp) .* ns.sfs;
        end
    else
        coss = ns.sfs .* cos(ns.phases);
        sins = sin(ns.phases) .* ns.sfs;
    end
    w = GetEval(Expt,'wi');
    p2d = atan(Expt.Stimvals.px/Expt.Stimvals.vd) * 180/pi;
    if isempty(xv)
        xv = -w/2:p2d:w/2;
    end
end


if isempty(Expt)
    if isempty(sfvals)
        sfvals = 1:size(ns.sfs,2);
    end
else
sf = GetEval(Expt,'sf');
df = Expt.Stimvals.f2 - sf;
nc = Expt.Stimvals.nc;
sfvals(1) = sf - df * nc/2;
sfvals(nc) = sfvals(1) + df * (nc-1);
if ismember(Expt.Stimvals.sM,[7 9 10])
  sfvals = exp(linspace(log(sfvals(1)),log(sfvals(nc)),nc));
else
  sfvals = linspace(sfvals(1),sfvals(nc),nc);
end
end
if isempty(ids)
    if verbose
        fprintf('Building id list');
    end
    ids = RlsRCc(Expt, ns.seeds, ns.ids, 'delays',delays, ...
		 'uselast','offset',seedoffset);
    res.rcid = ids;
end


if ismember(plottype,[1 2])
    subplot(1,1,1);
    for j = 1:length(delays)
        c = mean(coss(ids{j},:));
        s= mean(sins(ids{j},:));
        spc(j,:) = tospace(sfvals, xv, c,s);
        if isfield(ns,'dp')
            c = mean(rcoss(ids{j},:));
            s= mean(rsins(ids{j},:));
            rspc(j,:) = tospace(sfvals, xv, c,s);
        end
            
        ft(j,:) = abs(fft(spc(j,:)));
    end
    if plottype == 1
        imagesc(spc);
        title('xcorrelation');
    elseif plottype == 2
        imagesc(ft(:,1:20));
        title('FT of xcorr');
    end
    res.ft = ft;
    
elseif plottype == DPSF
    [nr,nc] = Nsubplots(length(delays));
    for j = 1:length(delays)
        for k = 1:length(ns.sfvals)
            sid = find(ns.sfs(ids{j},k));
            spc{j}(k,:) = mean(ns.dp(ids{j}(sid),:));
        end
        subplot(nr,nc,j);
        imagesc(spc{j});
    end
elseif plottype == MONOC
   
    for j = 1:length(delays)
        spc(j,:) = mean(ns.rsfs(ids{j},:));
        lspc(j,:) = mean(ns.lsfs(ids{j},:));
        if sfid
            id = find(ns.rsfs(ids{j},sfid) > 0);
            phases(j,:) = hist(ns.phases(ids{j}(id),sfid));
        end
    end
    
    if length(delays) > 1
        if plottype == 3
            imagesc(phases);
            res.phasehist = phases;
        else
            subplot(2,1,1);
            [X,Y,Z] = fillpmesh(sfvals,delays./10, spc);
            pcolor(X, Y, Z);
            shading('flat');
            title('R eye');
            
            subplot(2,1,2);
            [LX,LY,LZ] = fillpmesh(sfvals,delays./10, lspc);
            pcolor(LX, LY, LZ);
            shading('flat');
            title('L eye');
            xlabel('SF');
            ylabel('Delays');
            SetSubplots(2,1,[1:2],'caxis');
            if logscale(1)
                SetSubplots(2,1,[1:2],{'Xscale', 'Log', 'Xtick', [0.5 1 2 4 8]});
            end
            %imagesc(spc);
        end
    else
        plot(sfvals,spc);
    end
else
   
    for j = 1:length(delays)
        spc(j,:) = mean(ns.sfs(ids{j},:));
        if sfid
            id = find(ns.sfs(ids{j},sfid) > 0);
            phases(j,:) = hist(ns.phases(ids{j}(id),sfid));
        end
    end
    subplot(1,1,1);
    if length(delays) > 1
        if plottype == 3
            imagesc(phases);
            res.phasehist = phases;
        else
            [X,Y,Z] = fillpmesh(sfvals,delays, spc);
            pcolor(X, Y, Z);
            shading('flat');
            title('Binoc');
            %imagesc(spc);
        end
    else
        plot(sfvals,spc);
    end
end
res.means = spc;
if isfield(ns,'dp') & exist('rspc','var')
    res.rmeans = rspc;
end
res.delays = delays;
res.sfs = sfvals;
res.ids = ids;
res.sfvals = sfvals;

function sta = caltsta(im, ids)


function x = tospace(sfs, deg, c,s)


for j = 1:length(sfs)
   cmps(j,:) =  c(j) .* cos(deg * sfs(j) * 2 * pi) + s(j) .* sin(deg * sfs(j) * 2 * pi);
end
x = mean(cmps);