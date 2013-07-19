function [result, Expt] = PlotRevCorAny(Expt, varargin)
%
% result = PlotRevCorAny(Expt....)
% Plots up data for reverse correlation experiments. By default, a half Gaussian
% window with sd 2ms is used to smooth the sdfs, the variance betweeen sdfs
% as a function of time is estimated and time slices through the tuning
% curves are shown from the time when variance reaches 10% of peak to the
% time
% it falls back to 10% of peak
%
% PlotRevCorAny(Expt, 'sdfw', x) sets the Gaussian sd to x ticks (x/10 ms).
%
% PlotRevCorAny(Expt, 'slices', x) shows tuning curves for each of the
% delay times in vector x. 
%
% PlotRevCorAny(Expt, 'times', x) evaluates the sdfs for times x (relative to stimsulus transition)
% (default -100:10:1500). If sdfw is large, and new times are not given, a
% new default range is set
%
% For all of these arguments, times are in clock ticks, not ms 
%
% PlotRevCorAny(Expt, 'yval', Y) uses only stimuli for which the value of
%                                exp2 is in Y
% PlotRevCorAny(Expt, 'ytrial', Y) does the equivalent for subspace maps, seletcing
%                                  Trials for which exp2 is Y
%
% Second order plots come in several types:
% PlotRevCorAny(Expt, 'secondorder', [4 1]);   
%            is a standard second order plot, showing data for all
%            combinations
% PlotRevCorAny(Expt, 'secondorder', [4 n]);   
%            shows a reduced data set, showing data for all stimuli at t =
%            0,
%            for each stimulus at t = -1 whose index is n or greater.
%            this can be combined with a third element in the specification
%            vector to limit the number of comparions:
% PlotRevCorAny(Expt, 'secondorder', [4 n m]);
%            this shows comparisons for only m+1 values at t = -1. So
%
% PlotRevCorAny(Expt, 'secondorder', [4 3 0]); Shows only the responses for
% where the third stimulus preceded the stimulus at t = 0
%
% Other uses of the 'secondorder' vector
%  [5 n m] is much the same as 4, but selects a limited number of stimuli
%  at t=0, showing all precedints. So [4 3 0] shows all responses to the
%  third stimulus at t=0, with a separate line for each preceding stimulus
%  value.
%
%  [1 1]  for AC expts. Plot cases where preceding = AC and preceding =
%  Corr separately
%  [2 1] Another special case, plotting preceding disp < 0 and >=0
%  separately.
%  [3 1] Plot separate RCs for preceding x < current x and preceding x >
%  current x   [3 1 0] shows only precedingx < current, [3 2 0] shows only
%  preceding x > current.
%
% ...'collapse' ignores one dimension of a 2-D expt. 
%     ....'collapse',1 collapses over values of expt1, 2 collapses expt2
% Add uncorr to plot
% make plots for AC -> and COR ->
% get default slices from variance wrt time.
%
% For DCPD expt (CP for time varying disparity). use
% PlotRevCorAny(Expt,'ytrial',0,'psych'...

times = -100:10:1500;
np = 1;
sdfw = 20;
plotsum = 0;
plotbest = 0;
slices = [500 600 700 800 900 1000];
xrange = [];
labela = 'ACRCFiga';
labelb = 'ACRCFigb';
labelc = 'ACRCFigc';
makeall = 0;  %fill in sdf matrix regardless of count
secondorder = 0;
nloops = 1;
legendpos = 0;
result = [];
showplot = 1;
needvar = 0;
smtype = 'halfg';
autoslice = 1;
legendstyle = 0;
label_this_yval = 1;
interp = 0;
summarize = 0;
nmin = 200;
getpsych = 0;
npsych = 1;
profile = 0;
iskip = 1;
aloop = 1;
collapse = [0 0 0];
collapsetypes = {};
setyvals = [];
type = [];
btype = [];
ctype = [];
extraidx = {};
minplottime = 0;
autonmin = 0;
noextras = 0;
plotlfp = 0;
settimes = 0;
filltimes = 0;
framerate = 0;
addheader = 0;
autobox = 0;
nz = 0;
blankref = 0;

binsiz = [0 0];
centerRFvals = 0;
if profile
    tic;
end

showpcolor = 0;
findarg = {};
if ischar(Expt)  % Named a file, not passed an Expt structure
    file = Expt;
    load(file);
end

SpkDefs;

%
% Have to do calculations involvintg high signal trials first, before 
% low signal trials are selected out.
if isfield(Expt.Trials,'RespDir')
    if ~isfield(Expt.Trials,'Signal')
        if isfield(Expt.Trials,'ori')
            [Expt.Trials.Signal] = deal(Expt.Trials.ori);
        elseif isfield(Expt.Trials,'dx')
            [Expt.Trials.Signal] = deal(Expt.Trials.dx);
        end
    end
    if isfield(Expt.Trials,'Signal')
    ms = mean([Expt.Trials.Signal]);
        did = find([Expt.Trials.Dc] >0 & [Expt.Trials.Signal] > ms);
        uid = find([Expt.Trials.Dc] >0 & [Expt.Trials.Signal] < ms);
%CP is counts (RespDir >0) as pref.
%upstim is the stim associated with respdir > 0 (?actually a down saccade)
        result.cp.uprespdir = 1;
        if mean([Expt.Trials(uid).RespDir]) > mean([Expt.Trials(did).RespDir])
            result.cp.upstim = mean([Expt.Trials(uid).Signal]);
            result.cp.dnstim = mean([Expt.Trials(did).Signal]);
        else
            result.cp.upstim = mean([Expt.Trials(did).Signal]);
            result.cp.dnstim = mean([Expt.Trials(uid).Signal]);
        end
        highs = prctile([Expt.Trials.Dc],90);
        dsig = find([Expt.Trials.Dc] >= highs & [Expt.Trials.Signal] > ms);
        usig = find([Expt.Trials.Dc] >= highs & [Expt.Trials.Signal] < ms);
        if isfield(Expt.Trials,'LFP')
        result.cp.siglfpwr(1,:,:) = mean(abs(fft(cat(3,Expt.Trials(usig).LFP))),3);
        result.cp.siglfpwr(2,:,:) = mean(abs(fft(cat(3,Expt.Trials(dsig).LFP))),3);
        result.cp.siglfpt(1,:,:) = mean(cat(3,Expt.Trials(usig).LFP),3);
        result.cp.siglfpt(2,:,:) = mean(cat(3,Expt.Trials(dsig).LFP),3);
        end
    end
end

j = 1; 
while j < nargin
    if strncmpi(varargin{j},'addheader',4)
        addheader = 1;
    elseif strncmpi(varargin{j},'autoslices',4)
        autoslice = 1;
    elseif strncmpi(varargin{j},'bin',3)
        j = j+1;
        binsiz = varargin{j};
    elseif strncmpi(varargin{j},'blankref',7)
        blankref = 1;
    elseif strncmpi(varargin{j},'slices',4)
        j = j+1;
        slices = varargin{j};
        autoslice = 0;
    elseif strncmpi(varargin{j},'add',3)
        plotsum = 1;
    elseif strncmpi(varargin{j},'fbox',4)
        smtype = 'box';
        autobox = 1;
    elseif strncmpi(varargin{j},'box',3)
        smtype = 'box';
        minplottime = times(1) + sdfw;
    elseif strncmpi(varargin{j},'block',4)
        j = j+1;
        useblk = varargin{j};
        blk = Expt.Header.BlockStart;
        blk(end+1) = max([Expt.Trials.Trial]);
        trials = [];
        for k = 1:length(useblk)
            trials = [trials find(ismember([Expt.Trials.Trial],[blk(useblk(k)):blk(useblk(k)+1)]))];
        end
        Expt.Trials = Expt.Trials(trials);
    elseif strncmpi(varargin{j},'collapse',5)
        j = j+1;
        if ischar(varargin{j})
            collapsetypes{1} = varargin{j};
        elseif length(varargin{j}) > 1
        collapse = varargin{j};
        else
        collapse(varargin{j}) = 1;
        end
    elseif strncmpi(varargin{j},'range',3)
        j = j+1;
        xrange = varargin{j};
    elseif strncmpi(varargin{j},'makeall',7)
        makeall = 1;
    elseif strncmpi(varargin{j},'mint',3)
        j = j+1;
        minplottime = varargin{j};
    elseif strncmpi(varargin{j},'exp3',4)
        j = j+1;
        ctype = varargin{j};
    elseif strncmpi(varargin{j},'exp2',4)
        j = j+1;
        btype = varargin{j};
    elseif strncmpi(varargin{j},'exp',3)
        j = j+1;
        type = varargin{j};
    elseif strncmpi(varargin{j},'best',4)
        plotbest = 1;
    elseif strncmpi(varargin{j},'figa',4)
        j = j+1;
        labela = varargin{j};
    elseif strncmpi(varargin{j},'figb',4)
        j = j+1;
        labelb = varargin{j}; 
    elseif strncmpi(varargin{j},'figc',4)
        j = j+1;
        labelc = varargin{j}; 
    elseif strncmpi(varargin{j},'filltimes',5)
        filltimes = 1;
        j = j+1;
        framerate = varargin{j};
    elseif strmatch(varargin{j},'lfp')
      plotlfp = 1;
      if length(varargin) > j & isnumeric(varargin{j+1})
          plotlfp = varargin{j+1};
          j = j+1;
      end
    elseif strmatch(varargin{j},'legendpos')
      legendpos = varargin{j+1};
      j = j+1;
    elseif strncmpi(varargin{j},'nmin',4)
      nmin = varargin{j+1};
      if ischar(nmin) & strmatch(nmin,'auto')
          autonmin = 1;
      end
      j = j+1;
    elseif strncmpi(varargin{j},'noplot',5)
        showplot = 0;
    elseif strncmpi(varargin{j},'noextra',5)
        noextras = 1;
    elseif strncmpi(varargin{j},'pcolor',5)
        showpcolor = 2;
        if length(varargin) > j & isnumeric(varargin{j+1})
            showpcolor = varargin{j+1};
            j = j+1;
        end
    elseif strncmpi(varargin{j},'profile',5)
        profile = 1;
    elseif strncmpi(varargin{j},'psych',5)
        getpsych = 1;
    elseif strncmpi(varargin{j},'upchoice',3)
        getpsych = 2;
    elseif strncmpi(varargin{j},'downchoice',4)
        getpsych = 3;
    elseif strncmpi(varargin{j},'interp',5)
        interp = 1;
        if length(varargin) > j & isnumeric(varargin{j+1})
            interp = varargin{j+1};
            j = j+1;
        end
    elseif strncmpi(varargin{j},'sdfw',4)
        j = j+1;
        sdfw = varargin{j};
    elseif strncmpi(varargin{j},'summary',6)
        summarize = 1;
    elseif strncmpi(varargin{j},'secondorder',6)
        j = j+1;
        secondorder = varargin{j}(1);
        aloop = varargin{j}(2);
        if length(varargin{j}) > 2
            nloops = aloop+varargin{j}(3);
        elseif ismember(secondorder,[3 4 5])
            nloops = 0;
        elseif ismember(secondorder,[6])
            nloops = 1;
        elseif ismember(aloop,1)
            nloops = 2;
        end
    elseif strncmpi(varargin{j},'select',6)
        j = j+1;
        idx = eval(['find(' varargin{j} ');']);
        Expt.Trials = Expt.Trials(idx);
    elseif strncmpi(varargin{j},'skip',4)
        j = j+1;
        iskip = varargin{j};
    elseif strncmpi(varargin{j},'splitce',6)
        secondorder = 1;
        aloop = 1;
        nloops = 2;
    elseif strncmpi(varargin{j},'splitdx',6)
        secondorder = 2;
        aloop = 1;
        nloops = 2;
    elseif strncmpi(varargin{j},'times',5)
        j = j+1;
        times = varargin{j};
        settimes = 1;
    elseif strncmpi(varargin{j},'twoslice',5)
        if isempty(type)
            type = Expt.Stimvals.et;
        end
        if isempty(btype)
            btype = Expt.Stimvals.e2;
        end
        a = PlotRevCorAny(Expt,varargin{1:j-1},'exp',type,'exp2',btype,'collapse',2);
        b = PlotRevCorAny(Expt,varargin{1:j-1},'exp',btype,'exp2',type,'collapse',2);
        result = a;
        if size(b.sdfs.s,1) == size(a.sdfs.s,1)
            if isfield(a,'y')
                result.y = cat(2, a.y, b.y);
            end
            if isfield(a,'vars')
                result.vars = a.vars+b.vars;
            end
            result.yvals = cat(2, a.yvals, b.yvals);
            result.sdfs.s = cat(2,a.sdfs.s,b.sdfs.s);
            result.sdfs.x = cat(2,a.sdfs.x,b.sdfs.x);
            result.sdfs.n = cat(2,a.sdfs.n,b.sdfs.n);
            result.sdfs.y = cat(2,a.sdfs.y,b.sdfs.y);
            if isfield(a.sdfs,'lfp')
                result.sdfs.lfp = cat(2,a.sdfs.lfp,b.sdfs.lfp);
                result.sdfs.lfpn = cat(2,a.sdfs.lfpn,b.sdfs.lfpn);
            end
        else
            result.bresult = b;
        end
        GetFigure(labela); hold off;
        PlotRC(result,'label',labela);
        GetFigure(labelb); hold off;
        PlotRC(result,'sdfs','label',labelb);
        return;
    elseif strncmpi(varargin{j},'ymax',4)
        j = j+1;
        yv = Expt.Stimvals.et;
        for k = 1:length(Expt.Trials)
            yvals(k) = Expt.Trials(k).(yv)(end);
        end
        idx = find( yvals <= varargin{j} );
        Expt.Trials = Expt.Trials(idx);
    elseif strncmpi(varargin{j},'xtrial',4)
        j = j+1;
        yv = Expt.Stimvals.et;
        for k = 1:length(Expt.Trials)
            yvals(k) = Expt.Trials(k).(yv)(end);
        end
        idx = find(ismember(yvals,varargin{j}));
        Expt.Trials = Expt.Trials(idx);
    elseif strncmpi(varargin{j},'ytrial',4)
        j = j+1;
        yv = Expt.Stimvals.et;
        for k = 1:length(Expt.Trials)
            yvals(k) = Expt.Trials(k).(yv)(end);
        end
        idx = find(ismember(yvals,varargin{j}));
        Expt.Trials = Expt.Trials(idx);
    elseif strncmpi(varargin{j},'ztrial',4)
        j = j+1;
        yv = Expt.Stimvals.e3;
        for k = 1:length(Expt.Trials)
            yvals(k) = Expt.Trials(k).(yv)(end);
        end
        idx = find(ismember(yvals,varargin{j}));
        Expt.Trials = Expt.Trials(idx);
    elseif strncmpi(varargin{j},'trials',6)
        j = j+1;
        Expt.Trials = Expt.Trials(varargin{j});
    elseif strncmpi(varargin{j},'yval',4)
        j = j+1;
        setyvals = varargin{j};
    elseif strncmpi(varargin{j},'zval',4)
        nz = nz+1;
        zvarg(nz) = j;
        j = j+1;
        ztype{nz} = varargin{j};
        j = j+1;
        setzvals{nz} = varargin{j};
    end
    j = j+1;
end


if isfield(Expt.Trials,'Fr') & length(unique([Expt.Trials.Fr])) > 1
    Frs = unique([Expt.Trials.Fr]);
    for j = 1:length(Frs)
        id = find([Expt.Trials.Fr] == Frs(j));
        bExpt = Expt;
        bExpt.Trials = Expt.Trials(id);
        if j == length(Frs)
        subres{j} = PlotRevCorAny(bExpt,varargin{:});
        else
        subres{j} = PlotRevCorAny(bExpt,varargin{:},'figa',labelc);
        GetFigure(labelc);
        title(sprintf('Fr = %.0f',Frs(j)));
        end
        if isfield(subres{j},'bestdelay') && subres{j}.bestdelay <= size(subres{j}.y,3) %missing for LFP files.
            ys{j} = subres{j}.y(:,:,subres{j}.bestdelay);
        end
        subres{j}.Fr = Frs(j);
    end
    result = subres{j};
    result.subres = {subres{1:j-1}};
    result.ctype = 'Fr';
    result.name = Expt.Header.Name;
    GetFigure(labelb);
    hold off;
    PlotRC(result,'best','label',labelb);
    return;
end

if ~isempty(ctype)
    cs = unique([Expt.Trials.(ctype)]);
    if length(unique(cs)) > 1
    for j = 1:length(cs)
        id = find([Expt.Trials.(ctype)] == cs(j));
        bExpt = Expt;
        bExpt.Trials = Expt.Trials(id);
        if j == length(cs)
        subres{j} = PlotRevCorAny(bExpt,varargin{:});
        else
        subres{j} = PlotRevCorAny(bExpt,varargin{:},'figa',labelc);
        GetFigure(labelc);
        title(sprintf('%s = %.2f',ctype,cs(j)));
        end
        subres{j}.(ctype) = cs(j);
    end
    result = subres{j};
    result.subres = {subres{1:j-1}};
    result.name = Expt.Header.Name;
    return;
    end
end

if nz & size(setzvals{1},1) > 1 & size(setzvals{1},2) > 1
    for j = 1:size(setzvals{1},1)
        args = {};
        for k = 1:length(ztype)
            args = {args{:}, 'zval', ztype{k}, setzvals{k}(j,:)};
        end
        res{j} = PlotRevCorAny(Expt, varargin{1:zvarg(1)-1},args{:});
    end
    result = res{1};
    result.subres = res(2:end);
    return;
end


    if isfield(Expt.Header,'frameperiod') & ~isempty(Expt.Header.frameperiod) %number calculated from Spike2 events
        frameperiod = Expt.Header.frameperiod;
    elseif isfield(Expt.Stimvals,'fz') & ~isempty(Expt.Stimvals.fz)
        frameperiod = 10000/Expt.Stimvals.fz;
    else
        frameperiod = 104;
    end
    
if sdfw == 0
    if strcmp(smtype,'box')
        autobox = 1;
        sdwf = frameperiod;
    end
end

if sdfw > range(times)/5  && ~settimes
    times = -sdfw:sdfw/10:sdfw * 10;
end
result.calctime(1) = now;
colors = mycolors;
excolors = mycolors;
excolors{4} = [0.6 0.4 0]; %% dotted yellow invisible
colors = {colors{:} colors{:}}; % just in case second order makes many lines

if ~isfield(Expt.Header,'RCparams')
    if isfield(Expt.Header,'StoreErr')
        for j = 1:length(Expt.Trials)
            lens(j) = length(Expt.Trials(j).Start);
        end
        nf = max(lens);
        for j = 1:length(Expt.Trials)
            if(length(Expt.Trials(j).Start) < nf)
                Expt.Trials(j).Start(end:nf) = NaN;
            end
        end
        fprintf('NaN Paddding for Store Bug\n');
    else
        for j = 1:length(Expt.Trials)
            lens(j) = length(Expt.Trials(j).Start);
            id = find(Expt.Trials(j).Start < 0);
            if ~isempty(id)
                Expt.Trials(j).Start(id) = NaN;
            end
        end
        nf = min(lens);
        nf = floor(prctile(lens,10));
        if Expt.Stimvals.sM == 26  %bad fix stored too, so many short trials.
            nf = floor(prctile(lens,75));
        else
            nf = floor(prctile(lens,10));
        end
        Expt.Trials = Expt.Trials(find(lens >= nf));
        for j = 1:length(Expt.Trials)
            Expt.Trials(j).Start = Expt.Trials(j).Start(1:nf);
        end
    end
    

    starts = [Expt.Trials.Start];

    Expt.Header.Name = strrep(Expt.Header.Name,'\','/');
    result.name = Expt.Header.Name;
    if isfield(Expt.Header,'probe') || isfield(Expt,'probes')
        if isfield(Expt,'probes')
            result.probe = median(Expt.probes);
        else
        result.probe = Expt.Header.probe;
        end
        if isfield(Expt.Header,'Clusters')
            if length(Expt.Header.Clusters{1}) > 1
                result.Cluster = {Expt.Header.Clusters{1}{:,round(result.probe)}};
            else
                result.Cluster = Expt.Header.Clusters{1};
            end
        end
    else
        if isfield(Expt.Header,'Clusters')
            result.Cluster = Expt.Header.Clusters;
        end
    end



if isempty(type)
    if ischar(Expt.Stimvals.et)
        type = Expt.Stimvals.et;
    else
        if Expt.Stimvals.et == 263
            type = 'dO';
        elseif Expt.Stimvals.et == 237
            type = 'dO';
        end
    end
end
if isempty(btype)
    if ischar(Expt.Stimvals.e2)
        btype = Expt.Stimvals.e2;
    else
        if Expt.Stimvals.e2 == 112
            btype = 'ce';
        elseif Expt.Stimvals.e2 == 117
            btype = 'ce';
        end
    end
end
    if strcmp(btype,'Pd') & strcmp(type,'Dc')
        type = 'Pd';
        btype = 'Dc';
    end
    if strcmp(btype,'dx') & strcmp(type,'Dc')
        type = 'dx';
        btype = 'Dc';
    end
    if strcmp(btype,'or') & strcmp(type,'Dc')
        type = 'or';
        btype = 'Dc';
        [Expt.Trials.Signal] = deal(Expt.Trials.ori);
    end
    if isempty(btype) && ismember(Expt.Stimvals.sM,[13 14]) 
        if isfield(Expt.Trials,'ph')
            btype = 'ph';
        elseif isfield(Expt.Trials,'lph') & strcmp(type,'Ol')
            btype = 'lph';
        elseif isfield(Expt.Trials,'rph') & strcmp(type,'Or')
            btype = 'rph';
        end
    end
    if strcmp(btype,'e0') | collapse(2)
        if ~strcmp(Expt.Stimvals.e3,'e0') && isfield(Expt.Trials,Expt.Stimvals.e3) && collapse(3) == 0
            btype = Expt.Stimvals.e3;
        else
            btype = [];
        end
    end
%NB NOT types. So that it matches PlotExpt
    result.type{1} = type;
    result.type{2} = btype;
    if strcmp(btype,'ce')
        linestyles = {'--', '-', ':','-',':',};
        label_this_yval = 2;
    elseif getpsych == 1
        linestyles = {'--', '-', ':','-.','-','--',':','-.','-','--','--', '-', ':','-.','-','--',':','-.','-','--'};
        label_this_yval = 2;
    else
        linestyles = {'-', '--', ':','-.','-',':','--','-.','-','-', '--', ':','-.','-',':','--','-.','-'};
    end


    %
    % first make sure all trials have the same length vectors
    % describing the stimulus

    if isfield(Expt.Trials,'Spikes') & ~plotlfp
    res.spksum = 0;
    for j = 1:length(Expt.Trials)
        Expt.Trials(j).eval = Expt.Trials(j).(type)(end);
        len(j) = eval(['length(Expt.Trials(j).' type ');']);
        res.spksum = res.spksum + length(Expt.Trials(j).Spikes);
    end

    if res.spksum == 0
        if isfield(Expt.Header,'expname') & isfield(Expt.Header,'probe')
            fprintf('No Spikes in %s,%s,%d\n',Expt.Header.Name,Expt.Header.expname,Expt.Header.probe);
        else
            fprintf('No Spikes in %s\n',Expt.Header.Name);
        end
        return;
    end
    else
        for j = 1:length(Expt.Trials)
            Expt.Trials(j).eval = Expt.Trials(j).(type)(end);
            len(j) = eval(['length(Expt.Trials(j).' type ');']);
        end
    end
    if isfield(Expt.Header,'StoreErr')
        tlen = max(len);
        Nf = tlen;
    else
        tlen = floor(mean(len)+0.1);
% if nf is < tlen, need to use nf. Not sure why  check for diff < tlen/3        
        if nf < tlen & tlen -nf < tlen/3;
            tlen = nf;
        end
%what is stored in Trials.Nf is the number of frames completed _before_
%the last frame (becuase it is send then). So Add one.
%Spike times are all given relating to Trials.TrialStart. So referstnce
%stim frame times to this also (see FindTrials, which sets trigger).
        if isfield(Expt.Trials,'Nf')
            nfs = Expt.Trials.Nf;
            Nf = 1+round(mean(nfs(find(nfs > 0))));
        else
            if isfield(Expt.Stimvals,'Nf')
                Nf = Expt.Stimvals.Nf;
            else
                Nf = GetEval(Expt,'nf');
            end
        end
        if Nf < tlen & Nf > 0
            tlen = Nf;
        end
    end

    if profile
        fprintf('Filling trials at %.2f\n',toc);
    end
    for k = 1:length(collapsetypes)
        collapseval(k) = prctile(cat(1, Expt.Trials.(collapsetypes{k})),50);
    end


    if autobox
        sdfw = frameperiod;
        minplottime = times(1) + sdfw;
    end
    lastf = tlen;
    for j = 1:length(Expt.Trials)
        if ~isempty(btype) & length(Expt.Trials(j).(btype)) == 1
            Expt.Trials(j).(btype)(1:tlen,1) = Expt.Trials(j).(btype);
        end
        if len(j) < tlen
            Expt.Trials(j).(type)((len(j)+1):tlen) = NaN;
            if isfield(Expt.Trials,'ce')
                Expt.Trials(j).ce((len(j)+1):tlen) = NaN;
            end
            if ~isempty(btype)
                Expt.Trials(j).(btype)((len(j)+1):tlen) = NaN;
            end
            if isfield(Expt.Trials,'st')
                Expt.Trials(j).st((len(j)+1):tlen) = NaN;;
            end
            if isfield(Expt.Trials,'me')
                Expt.Trials(j).me((len(j)+1):tlen) = NaN;;
            end
            lastf = len(j);
        elseif len(j) > tlen
            Expt.Trials(j).(type) = Expt.Trials(j).(type)(1:tlen);
            if ~isempty(btype)
                if length(Expt.Trials(j).(btype)) >= tlen
                Expt.Trials(j).(btype) = Expt.Trials(j).(btype)(1:tlen);
                else
                    l = length(Expt.Trials(j).(btype));
                    Expt.Trials(j).(btype)(l+1:tlen) = NaN;
                end
            end
            if isfield(Expt.Trials,'ce') & length(Expt.Trials(j).ce) > 1
                Expt.Trials(j).ce  = Expt.Trials(j).ce(1:tlen);
            end
            if isfield(Expt.Trials,'st')
                Expt.Trials(j).st = Expt.Trials(j).st(1:tlen);
            end
            if isfield(Expt.Trials,'me')
                if length(Expt.Trials(j).me) >= tlen
                Expt.Trials(j).me = Expt.Trials(j).me(1:tlen);
                else
                Expt.Trials(j).me = ones(1,tlen) .* Expt.Trials(j).me;
                end
            end
        end
 % sometimes
        if strcmp(type,'or')
            id = find(Expt.Trials(j).or > -1000);
           Expt.Trials(j).or(id) = mod(Expt.Trials(j).or(id),180);
        end
        id = find(isnan(Expt.Trials(j).Start));
        id = id(find(id <= tlen));
        for k = 1:length(id)
            Expt.Trials(j).(type)(id(k)) = Expt.Trials(j).(type)(id(k)-1);
        end
        if filltimes
            Expt.Trials(j).Start = Expt.Trials(j).Start(1) + framerate .* [0:tlen-1]'; 
            Expt.Trials(j).End = Expt.Trials(j).Start + framerate; 
        end
        for k = 1:length(collapsetypes)
            Expt.Trials(j).(collapsetypes{k})(1:tlen) = collapseval(k);
        end
        id = find((Expt.Trials(j).(type)) > 1000);
        Expt.Trials(j).(type)(id) = -10000;
        tstart = Expt.Trials(j).Start(1) - Expt.Trials(j).TrialStart;
        tend = Expt.Trials(j).Start(lastf) - Expt.Trials(j).TrialStart;
        if isfield(Expt.Trials,'Spikes');
        Expt.Trials(j).Count = sum(Expt.Trials(j).Spikes > tstart & Expt.Trials(j).Spikes < tend);
        end
        if isfield(Expt.Trials,'Dc')
        Expt.Trials(j).Dc = Expt.Trials(j).Dc(end);
        end
        Expt.Trials(j).TrialDur = tend-tstart;
    end
    rcparams.btype = btype;
    rcparams.type = type;
else
    btype = Expt.RCparams.btype;
    type = Expt.RCparams.type;
    starts = [Expt.Trials.Start];
end

if(autoslice)
    needvar = 1;
end


if profile
    fprintf('Done at %.2f\n',toc);
end


xvals = unique([Expt.Trials.(type)]);
if isfield(Expt.Trials,'st') & secondorder & ~ismember(-1009,xvals)
    for j = 1:length(Expt.Trials)
        id = find([Expt.Trials(j).st] == 0);
        Expt.Trials(j).(type)(id) = -1009;
    end
    xvals = unique([Expt.Trials.(type)]);
end
xvals = xvals(find(~isnan(xvals)));
if ~isempty(btype)
    yid = find(~isnan([Expt.Trials.(btype)]));
    ymat = [Expt.Trials.(btype)];
    %    yvals = unique([Expt.Trials.(btype)](yid));
    yvals = unique(ymat(yid));
else
    yvals = [];
end

if binsiz(1)
    newx = min(xvals)+binsiz(1)/2:binsiz(1):max(xvals);
    for j = 1:length(xvals)
        [a,b] = min(abs(xvals(j)-newx));
        newid(j) = b;
    end
    for j = 1:length(Expt.Trials)
        [a, Expt.Trials(j).xid] = ismember(Expt.Trials(j).(type), xvals);
        Expt.Trials(j).(type) = newx(newid(Expt.Trials(j).xid));
    end
    xvals = newx;
end

if length(binsiz) > 1 & binsiz(2)
    newy = min(yvals)+binsiz(2)/2:binsiz(2):max(yvals);
    for j = 1:length(yvals)
        [a,b] = min(abs(yvals(j)-newy));
        newid(j) = b;
    end
    for j = 1:length(Expt.Trials)
        [a, Expt.Trials(j).yid] = ismember(Expt.Trials(j).(btype), yvals);
        Expt.Trials(j).(btype) = newx(newid(Expt.Trials(j).yid));
    end
    yvals = newy;
end

if ~isempty(setyvals)
    if sum(ismember(yvals, setyvals)) == 0 %% no  matches
        [a,b] = min(abs(yvals - setyvals(1)));
        yvals = yvals(b);
    else
        yvals = setyvals;
    end
end

result.yvals = yvals;
if ~isempty(xrange)
    idx = find(xvals >= xrange(1) & xvals <= xrange(2));
    xvals = xvals(idx);
end

if nloops == 0
    if ismember(secondorder, [4 5])
        nloops = length(xvals);
    elseif secondorder == 3
        nloops = aloop+1;
    else
        nloops = 0;
    end
end


endoffset = 5; %% 0.5ms

if isfield(Expt,'RCparams')
    tx = [Expt.Trials.(type)]; %safe if RCparams - means this is checked already
    ty = [Expt.Trials.(btype)];
    durs = [Expt.Trials.dur];
else
    for j = 1:length(Expt.Trials)
        tx(:,j) = Expt.Trials(j).(type);
        if ~isempty(yvals)
            y = Expt.Trials(j).(btype);
            if length(y) < nf-1
            ty(:,j) = NaN;
            ty(1:length(y),j) = Expt.Trials(j).(btype);
            else
            ty(:,j) = Expt.Trials(j).(btype);
            end
        end
        durs(:,j) = [diff(Expt.Trials(j).Start); endoffset + Expt.Trials(j).End(end) - Expt.Trials(j).Start(end)];
        Expt.Trials(j).dur = durs(1:tlen,j);
    end
end

xv = xvals(find(xvals > -1000));
xbin = {};
k = 1;
if min(diff(xv)) < range(xv)/50 & max(diff(xv)) > range(xv)/25 %% rounding issues
    id = find(diff(xvals) < range(xv)/50);
    x = xvals(id(1));
    for j = 1:length(id)
        xi = find(tx == xvals(id(j)+1));
        if j > 1 && diff(id(j-1:j)) > 1
            x = xvals(id(j));
            xbin{k} = tx(xi);
            k = k+1;
        end
        tx(xi) = x;
    end
    xvals = unique(tx);
end
result.xbinned = xbin;
%
%
% check Stimulus durations and make sure that they are all the same.
% Exclude any outliers (more that 30% of frame period away from mean)
%

if isfield(Expt,'RCparams')
    stimdur = Expt.RCparams.stimdur;
else
    if profile
        fprintf('getting duriation at %.2f\n',toc);
    end
    stimdur = mode(durs(:));
    rcparams.stimdur = stimdur;
    if profile
        fprintf('Done at %.2f\n',toc);
    end
end

result.nbad = length(find(abs(durs(2:end,:) - stimdur) > frameperiod/3));
result.nframes = prod(size(tx));
result.frameperiod = frameperiod;
result.stimdur = stimdur;
result.Nf = Nf;
result.Stimvals.et = type;
result.Stimvals.e2 = btype;
result.sdfw = sdfw;
if ~isempty(ctype)
    result.ctype = ctype;
end
result.Header = Expt.Header;

badtimes = find(abs(durs - stimdur) > frameperiod/3 | isnan(durs));
good(badtimes) = 0;
good(1,:) = 1;
good = ones(size([Expt.Trials.dur]));
if(iskip)
    good(1:iskip,:) = 0;
end

for z = 1:nz
    result.ztypes = ztype;
    result.zvals = setzvals;
    for j = length(Expt.Trials):-1:1
        tz = Expt.Trials(j).(ztype{z});
        id = find(tz < setzvals{z}(1) | tz > setzvals{z}(2));
        good(id,j) = 0;
        tzall(:,j) = tz;
    end
    if ~isfield(result,'ctype')
        result.ctype = ztype{1};
        result.(ztype{1}) = mean(tzall(find(good)));
    end
end

if autonmin
    txvals = unique(tx(:));
    tyvals = unique(ty(:));
    nstim = length(txvals) * length(tyvals);
    nmin = sum(good(:))./(nstim * 3);
    if getpsych
        nmin = nmin/3;
    end
end
result.nmin = nmin;

if profile
    fprintf('Start Respdir at %.2f\n',toc);
end

if getpsych & isfield(Expt.Trials,'RespDir')
    respdirs = [Expt.Trials.RespDir];
    id = find(respdirs == -1);
    %    good(:,id) = 1;
    for j = 1:length(respdirs)
        if respdirs(j) == 1
            good(find(good(:,j)),j) = 2;% Down Choices
        elseif respdirs(j) == -1
            good(find(good(:,j)),j) = 1; %Up Choices
        elseif respdirs(j) == 0
            good(:,j) = 0;
        end
    end
    if isfield(Expt.Trials,'Dc') 
        ms = mean([Expt.Trials.Signal]);
        dcvals = unique([Expt.Trials.Dc]);
        uid = find([Expt.Trials.Dc] == 0 & [Expt.Trials.RespDir] > 0);
        did = find([Expt.Trials.Dc] == 0 & [Expt.Trials.RespDir] < 0);
        if isfield(Expt.Trials,'Count')
        gcpnmin = 10;
        result.cp.cp = CalcCP([Expt.Trials(uid).Count],[Expt.Trials(did).Count]);
        result.cp.n = [length(uid) length(did)];
        a = mean([Expt.Trials([uid did]).Count]);
        b = std([Expt.Trials([uid did]).Count]);
        zup = ([Expt.Trials(uid).Count]-a)/b;
        zdn = ([Expt.Trials(did).Count]-a)/b;
        gcn = [];
        gcv = [];
        did = find([Expt.Trials.Dc] >0 & [Expt.Trials.Signal] > ms);
        uid = find([Expt.Trials.Dc] >0 & [Expt.Trials.Signal] < ms);
%CP is counts (dir >0) as pref.    
        for j = 2:length(dcvals)
            uid = find([Expt.Trials.Dc] == dcvals(j) & [Expt.Trials.RespDir] > 0 & [Expt.Trials.Signal] > ms);
            did = find([Expt.Trials.Dc] == dcvals(j) & [Expt.Trials.RespDir] < 0 & [Expt.Trials.Signal] > ms);
            if length(did) > gcpnmin & length(uid) > gcpnmin && abs(length(did)/(length(did) + length(uid)) -0.5) < 0.8
                a = mean([Expt.Trials([uid did]).Count]);
                b = std([Expt.Trials([uid did]).Count]);
                zup = [zup ([Expt.Trials(uid).Count]-a)/b];
                zdn = [zdn ([Expt.Trials(did).Count]-a)/b];
                gcv = [gcv dcvals(j)];
                gcn = [gcn [length(uid); length(did)]];
            end
            uid = find([Expt.Trials.Dc] == dcvals(j) & [Expt.Trials.RespDir] > 0 & [Expt.Trials.Signal] < ms);
            did = find([Expt.Trials.Dc] == dcvals(j) & [Expt.Trials.RespDir] < 0 & [Expt.Trials.Signal] < ms);
            if length(did) > gcpnmin & length(uid) > gcpnmin && abs(length(did)/(length(did) + length(uid)) -0.5) < 0.8
                a = mean([Expt.Trials([uid did]).Count]);
                b = std([Expt.Trials([uid did]).Count]);
                zup = [zup ([Expt.Trials(uid).Count]-a)/b];
                zdn = [zdn ([Expt.Trials(did).Count]-a)/b];
                gcv = [gcv dcvals(j) * -1];
                gcn = [gcn [length(uid); length(did)]];
            end
        end
        result.cp.gcp = CalcCP(zup,zdn);
        result.cp.gcpn = [length(zup) length(zdn)];
        result.cp.gvals = gcv;
        result.cp.gvaln = gcn;
        elseif plotlfp
           
            if isfield(Expt.Trials,'FTlfp')
                result.cp.lfpwr(1,:,:) = mean(cat(3,Expt.Trials(uid).FTlfp),3);
                result.cp.lfpwr(2,:,:) = mean(cat(3,Expt.Trials(did).FTlfp),3);
                len = size(result.cp.lfpwr,2);
                result.ftfrq = (0:len-1)/(len * Expt.Header.LFPsamplerate);
            else
            result.cp.lfpwr(1,:,:) = mean(abs(fft(cat(3,Expt.Trials(uid).LFP))),3);
            result.cp.lfpwr(2,:,:) = mean(abs(fft(cat(3,Expt.Trials(did).LFP))),3);
            len = size(result.cp.lfpwr,2);
            result.ftfrq = (0:len-1)/(len * Expt.Header.LFPsamplerate);
            end
            result.cp.lfpt(1,:,:) = mean(cat(3,Expt.Trials(uid).LFP),3);
            result.cp.lfpt(2,:,:) = mean(cat(3,Expt.Trials(did).LFP),3);

        end
        
    end
end

if summarize
    result.xv = unique(tx(find(~isnan(tx))));
    for j = 1:length(result.xv)
        result.nx(j) = length(find(tx == result.xv(j)));
    end
    if exist('ty')
        result.yv = unique(ty(find(~isnan(ty))));
    end
end


if showplot
    result.figa = GetFigure(labela);
    result.labela = labela;
    result.labelb = labelb;
    hold off;
end
h = [];

if size(xvals,1) > size(xvals,2)
    xvals = xvals';
end


if isempty(yvals)
    yvals = 0;
    ty = zeros(size(tx));
end

sdfs.extras = {};
sdfs.extraval = [];
isextra = zeros(size(tx));
nextra = 0;
nstims = prod(size(tx));
if isfield(Expt.Trials,'ce') & size([Expt.Trials.ce],1) > 1
    idx = find([Expt.Trials.ce] == 0);
    if ~isempty(idx) & length(idx) < nstims/length(yvals)
        nextra = nextra+1;
        extraidx{nextra} = idx;
        extralabel{nextra} = 'Uncorr';
        extraval(nextra) = IUNCORR;

    else
        uidx = idx;
    end
end

if isfield(Expt.Trials,'st')
    idx = find([Expt.Trials.st] == 0);
    if ~isempty(idx)
        nextra = nextra+1;
        extraidx{nextra} = idx;
        extralabel{nextra} = 'Blank';
        extraval(nextra) = IBLANK;
    end
end

if isfield(Expt.Trials,'me') && isempty(strmatch('me',{Expt.Stimvals.et Expt.Stimvals.e2 Expt.Stimvals.e3}))
    idx = find([Expt.Trials.me] == -1);
    if ~isempty(idx) && length(idx) < prod(size([Expt.Trials.me]))/3;
        nextra = nextra+1;
        extraidx{nextra} = idx;
        extralabel{nextra} = 'Left';
        extraval(nextra) = LMONOC;
    end
    idx = find([Expt.Trials.me] == 1);
    if ~isempty(idx) & length(idx) < prod(size([Expt.Trials.me]))/3;
        nextra = nextra+1;
        extraidx{nextra} = idx;
        extralabel{nextra} = 'Right';
        extraval(nextra) = RMONOC;
    end
end

idx = find(tx == -1006);
if ~isempty(idx)
    nextra = nextra+1;
    extraidx{nextra} = idx;
    extralabel{nextra} = 'Broad';
    extraval(nextra) = IBROADBAND;
end



if noextras
    nextra = 0;
    extraidx = {};
else
    if makeall
        xvals = setdiff(xvals,extraval);
        xvals = xvals(xvals > -1000);
    end
end
allextras = [];
for j = 1:nextra
    extraidx{j} = reshape(extraidx{j},length(extraidx{j}),1);
    allextras = [allextras; extraidx{j}];
    isextra(extraidx{j}) = j;
end

tmpy = ty;
tmpy(allextras) = NaN;
yvals = unique(tmpy(find(~isnan(tmpy))));
if isempty(yvals) %% all seem to be extras
    yvals = unique(ty);
end
if ~isempty(setyvals)
    if sum(ismember(yvals, setyvals)) == 0 %% no  matches
        [a,b] = min(abs(yvals - setyvals(1)));
        yvals = yvals(b);
    else
        yvals = setyvals;
    end
end

if ~isfield(Expt,'RCparams')
    Expt.RCparams = rcparams;
end

loopctr = 1;

if profile
    fprintf('Structure ready at %.2f\n',toc);
end


%
% good((i,j) == 2 means RespDir == 1, means Downward eye movement.
% good((i,j) == 1 means RespDir == -1, means upward eye movement.
%
% respir -1 is evaluated first

if size(good,1) > size(tx,1)
    good = good(1:size(tx,1),:);
end
ntg = 1;
if getpsych == 1
    respdirs = [-1 1];
    npsych = 2;
elseif getpsych == 2
    findarg = {'RespDir' 1};
    npsych = 1;   
elseif getpsych == 3
    findarg = {'RespDir' 2};
    respdirs = [1 -1];
    npsych = 1;  
end
pid = find(times > minplottime);

for j = 1:length(Expt.Trials)
    alltriggers{j} = [];
    if plotlfp
        nchs(j) = size(Expt.Trials(j).LFP,2);
    end
end

if plotlfp
    nch = max(nchs);
    Expt.Trials = Expt.Trials(find(nchs == nch));
end
for loop = aloop: nloops;
    nx = 1;
    if secondorder
        findarg{1} = 'secondorder';
        findarg{2}(1) = secondorder;
        findarg{2}(2) = loop;
        if ismember(secondorder,[4 5])
            findarg{2}(3) = xvals(loop);
            sxv = xvals(loop);
        end
    end
    for x = xvals;
        ny = 1;
        if secondorder == 6  % just get consec presnetations
            findarg{2}(3) = x;
            sxv = x;
        end
            
        for iy = 1:length(yvals) * npsych %% not fussy about row/column vector

            y = yvals(mod(iy-1,length(yvals))+1);
            if npsych > 1
                choiceval = ceil(iy/length(yvals));
                [Expt, tidx, n] = FindTrials(Expt, x, y, tx, ty, good, isextra, findarg{:},'RespDir',choiceval);
            else
                [Expt, tidx, n] = FindTrials(Expt, x, y, tx, ty, good, isextra, findarg{:});
                choiceval = 1;
            end
            %tidx is a list of Trials, containing at least one of the desired stim
            %types. n is the total number of stim presentations.
            for k =1:length(tidx)
                    alltriggers{tidx(k)} = [ alltriggers{tidx(k)} Expt.Trials(tidx(k)).Trigger];
            end
            if ~isempty(tidx) & (n > nmin);
                if plotlfp
                    [sdfs.lfp{nx,ny,loopctr}, a] = TrigLFP(Expt.Trials(tidx),[-500 2000],Expt.Header.LFPsamplerate,nch);
                    sdfs.lfptimes = a.lfptimes;
                    sdfs.lfpn(nx,ny,loopctr) = mean(a.n);
                    if isfield(Expt.Trials,'Spikes');
                        [sdf, n] = trigsdfa(Expt.Trials(tidx),sdfw,times,smtype);
                    else
                        sdf = 0;
                    end
               else
                    [sdf, n] = trigsdfa(Expt.Trials(tidx),sdfw,times,smtype);
                end
               
      
                if npsych > 1
                    sdfs.respdir(ny) = respdirs(choiceval);
                end
                sdfs.x(nx, ny, loopctr) = x;
                sdfs.y(nx, ny, loopctr) = y;
                sdfs.n(nx, ny, loopctr) = n;
                if secondorder
                    sdfs.z(nx, ny, loopctr) = sxv;
                end
                sdfs.s{nx, ny, loopctr} = sdf;
                if ~showplot
                elseif ~isfield(Expt.Trials,'Spikes')
                        hn = plot(sdfs.lfptimes./10,sdfs.lfp{nx,ny,loopctr}(:,plotlfp),'color',colors{nx},'linestyle',linestyles{mod(ny-1,4)+1});
                        hold on;
                elseif ny > 1
                    if( ~plotsum)
                        hn = plot(times(pid)/10,sdf(pid),'color',colors{nx},'linestyle',linestyles{mod(ny-1,4)+1});
                        hold on;
                    end
                elseif ny == 1
                    if ~plotsum
                        hn = plot(times(pid)/10,sdf(pid),'color',colors{nx},'linestyle',linestyles{loopctr});
                        hold on;
                    end
                elseif  ~isempty(tidx)
                    if ~plotsum
                        hn = plot(times(pid)/10,sdf(pid),'--','color','k');
                        hold on;
                    end
                end
                if (ny == label_this_yval | legendstyle ==1) & showplot
                    h(np) = hn;
                    if secondorder
                        labels{np} = sprintf('%.2f->%.2f n = %d',sxv,x,n);
                    else
                        if choiceval > 1
                            labels{np} = sprintf('%.2f n = %d,%d',x,n,sdfs.n(np));
                        else
                            labels{np} = sprintf('%.2f n = %d',x,n);
                        end
                    end
                    np = np+1;
                end
                ny = ny+1;
            else %not enough data
                if ny == 2 || (makeall && x > -1000)
                    sdfs.n(nx,ny) = n;
                    nx = nx+1;
                end
            end
        end
        %ny only incremented if data found. ny == 1 means no data found.
        if ny > 1
            nx = nx+1;
        end
    end
    loopctr = loopctr+1;
end



if profile
    fprintf('Main loop done at %.2f\n',toc);
end
nex = 1;
goodextras = 0;
for loop = aloop: nloops;
for j = 1:nextra * npsych
    if secondorder
        findarg{1} = 'secondorder';
        findarg{2}(1) = secondorder;
        findarg{2}(2) = loop;
        if ismember(secondorder,[4 5])
            findarg{2}(3) = xvals(loop);
        end
    end

    exarg = 'extra';
    if length(setyvals)
                exarg = 'extrabyy';
    end
    if npsych > 1
        lid = 1+ mod(j-1,npsych);
        choiceval = ceil(j/nextra);
        [Expt, tidx, n] = FindTrials(Expt, x, y, tx, ty, good, isextra, exarg, lid, findarg{:},'RespDir',choiceval);
    else
        [Expt, tidx, n] = FindTrials(Expt, x, y, tx, ty, good, isextra, exarg, j, findarg{:});
        lid = j;
        choiceval = loop;
    end
    for k =1:length(tidx)
        alltriggers{tidx(k)} = [ alltriggers{tidx(k)} Expt.Trials(tidx(k)).Trigger];
    end
    if ~isempty(tidx) & ((n > nmin) || makeall)
        if choiceval == 1
            goodextras = goodextras+1;
        end
        
        if plotlfp
            [sdfs.extras{lid,choiceval}.lfp, a] = TrigLFP(Expt.Trials(tidx),[-500 2000],Expt.Header.LFPsamplerate, nch );
            sdfs.extras{lid,choiceval}.lfpn = mean(a.n);
            sdfs.extras{lid,choiceval}.label = extralabel{lid};
        end
                
        if isfield(Expt.Trials,'Spikes')        
        [sdf, n] = trigsdfa(Expt.Trials(tidx),sdfw,times,smtype);
        if showplot
            h(np+nex-1) = plot(times(pid)/10,sdf(pid),':','color',excolors{nex - (choiceval-1)*goodextras},'linestyle',linestyles{3+choiceval-1},'linew',2);
            labels{np+nex-1} = [extralabel{lid} sprintf(' n = %d',n)];
        end
        sdfs.extras{lid,choiceval}.sdf = sdf;
        sdfs.extras{lid,choiceval}.n = n;
        sdfs.extraval(lid,choiceval) = extraval(lid);
        sdfs.extras{lid,choiceval}.label = extralabel{lid};
        elseif showplot
            h(np+nex-1) = plot(sdfs.lfptimes/10,sdfs.extras{lid,choiceval}.lfp(:,plotlfp),':','color',excolors{nex - (choiceval-1)*goodextras},'linestyle',linestyles{3+choiceval-1},'linew',2);
            labels{np+nex-1} = [extralabel{lid} sprintf(' n = %d',n)];
        end
        nex = nex+1;
    end
end
np = np+nex;
end

while length(sdfs.extras) > 0 && isempty(sdfs.extras{1})
    sdfs.extras = {sdfs.extras{2:end}};
end

if secondorder
for j = 1:length(alltriggers)
    alltriggers{j} = unique(alltriggers{j});
end
end

result.alltriggers = alltriggers;
if plotlfp
    [sdfs.alllfp] = TrigLFP(Expt.Trials,[-500 2000],Expt.Header.LFPsamplerate, nch,'triggers',alltriggers);    
    [sdfs.triallfp] = TrigLFP(Expt.Trials,[-500 2000],Expt.Header.LFPsamplerate, nch,'triggers',zeros(length(Expt.Trials),1));    
end
if profile
    fprintf('Extras done at %.2f\n',toc);
end

np = 1;

bestvar = 0;

sampleid = [];

if ~autoslice %calc var for slices
    needvar = 1;
end
if (plotbest  | needvar) & isfield(sdfs,'s')
    nv = 1;
    id = find(times < 0);
    if isempty(id)
        id(1) = 1;
    end
    bestnv = NaN;
    if ~autoslice
        for j = 1:length(slices)
            [diffs, id] = min(abs(slices(j) - times));
            ts(j) = id;
        end
    else
        ts = id(end):2:length(sdfs.s{1,1});
    end
    for j = ts;
        x = [];
        y = [];
        z = [];
        for k = 1:size(sdfs.s,1)
            x(k) = sdfs.x(k,1);
            for co = 1:size(sdfs.s,2)
                y(k,co) = 0; %may be summed across ex 3 below
                nc = size(sdfs.s,3);
                if nc == 0 %empty sdf = not enough trials
                    y(k,co) = NaN; 
                end
% for second order stimuli, just use first order kernels for variance
% estimate, so collapse across dimension 3 of sdfs
                for zi = 1:nc
                    if length(sdfs.s{k,co,zi}) >= j
                        y(k,co) = y(k,co) + sdfs.s{k,co,zi}(j)/nc;
                    elseif nc == 1
                        y(k,co) = NaN;
                    end
                end
            end
        end
        dvar = sum(var(y(~isnan(y))));
        vars(nv) = dvar;
        delays(nv) = times(j);
        sampleid(nv) = j;
        result.x(:,nv) = x;
        for co = 1:size(sdfs.s,2)
            result.y(:,:,nv) = y;
        end
        if(dvar > bestvar)
            bestvar = dvar;
            bestj = j;
            bestx = x;
            besty = y;
            bestnv = nv;
        end
        nv = nv + 1;
    end
    if exist('delays','var') %%missing if just doing lfp
    result.delays = delays;
    result.bestdelay = bestnv;
    result.vars = vars;
    else
        plotbest = 0;
    end
    result.timeoff = times(1);
    if showplot & plotbest & bestnv > 0
        GetFigure(labelb);
        hold off;
        for co = 1:size(sdfs.s,2)
            result.h(co) = plot(bestx,besty(:,co),'-','color',colors{np});
            hold on;
            np = np+1;
        end
        title(sprintf('%s at %d',Expt.Header.Name,bestj));
    end
    if showplot & exist('vars','var')
        GetFigure(labela);
        plot(delays./10,sqrt(vars).*4,'k','linewidth',2);
    end
else
    bestnv = 1;
end

if length(sampleid)
w = 100/((times(2)-times(1))* mean(diff(sampleid))); %10ms
else
w = 100/((times(2)-times(1))* 2); %10ms
end
if w < 1
    w = 1;
end
if exist('vars') & w > length(vars)/2;
    w = floor(length(vars)/2);
end
if isnan(bestnv) %LFP
    result.varratio(1) = 0;
elseif ~exist('vars','var')
    result.varratio = [1 0 0];
elseif bestnv <= w | bestnv >= length(vars)-w || isnan(w) 
    result.varratio(1) = 1;
    result.varratio(2) = 0;
    result.varratio(3) = 0;
elseif bestnv+w > length(vars)
    result.varratio(2) = mean(vars(end-2*w:end));
    result.varratio(3) = mean(vars(1:w*2));
    result.varratio(1) = result.varratio(2)/result.varratio(3);
elseif isnan(bestnv)
    result.varratio(1:3) = NaN;
else
    result.varratio(2) = mean(vars(bestnv-w:bestnv+w));
    result.varratio(3) = mean(vars(1:w*2));
    result.varratio(1) = result.varratio(2)/result.varratio(3);
end
if bestnv > w
    res.varratio(4) = prctile(vars,30);
end

if showpcolor
    hold off;
    startt = find(times > sdfw/2);
    if blankref && sum(sdfs.extraval == IBLANK)
        ib = find(sdfs.extraval == IBLANK);
    PlotPcolor(sdfs,times/10,[],interp,startt(1),'refsdf',sdfs.extras{ib}.sdf);
    else
    PlotPcolor(sdfs,times/10,[],interp,startt(1));
    end
    h = [];
end

if showplot
    if(legendpos < 6 & ~isempty(h))
        hid = find(ishandle(h) & h > 0);
        result.legend = legend(h(hid),labels{hid},legendpos);
    end
    title(sprintf('%s V%.1f',Expt.Header.Name,result.varratio(1)));
    ylabel('Rate (sp/s)');
    xlabel('Time (ms)');

    GetFigure(labelb);
    hold off;
end
h = [];
labels = {};


if ~plotbest & isfield(sdfs,'s') & bestnv > 0
    if showplot
        result.figb = GetFigure(labelb);
    end
    if autoslice
        [mv, mj] = max(result.vars);
        minv = min(result.vars);
        th = minv + (mv-minv)/10;
        id = find(result.vars > th);
        step = range(result.delays(id))/6;
        slices = round(result.delays(id(1)):step:result.delays(id(end)));
        slices = round(result.delays(mj)-2*step:step:result.delays(mj)+2*step);
        result.slices = slices;
        step = range(id)/6;
        slices = round(mj-step*2:step:mj+step*2);
        %cslices are coarse time samples, at the peak and at half height.
        id = find(result.vars > minv +(mv-minv)/2);
        cslices(2) = result.bestdelay;
        cslices(1) = id(1);
        cslices(3) = id(end);
    else
        %
        %
        result.slices = slices;
        for j = 1:length(slices)
            [diffs, id] = min(abs(slices(j) - times));
            sliceid(j) = id;
            sampleid(id) = id;
            delays(id) = times(id);
        end
        slices = sliceid;
        cslices = sliceid;
    end

    if profile
        fprintf('Slices ready at %.2f\n',toc);
    end

    %
    % id() keeps track of the data sample points where var was evaluated.
    % so sampleid(id(k)) is the data sample in the sdf associated with vars(k),
    % at time delays(id(k))
    %
    %
    if showplot
        subplot(1,1,1);
        hold off;
        for j = cslices;
            x = [];
            y = [];
            z = [];
            sn = sampleid(j);
            for k = 1:size(sdfs.s,1)
                if sdfs.x(k,1) > -1000
                x(k) = sdfs.x(k,1);
                for co = 1:size(sdfs.s,2)
                    if isempty(sdfs.s{k,co})
                        y(k,co) = NaN;
                    else
                        y(k,co) = sdfs.s{k,co}(sn);
                    end
                end
                else
                    x(k) = NaN;
                    y(k,co) = NaN;
                end
            end
            for co = 1:size(sdfs.s,2)
                ls = 1+mod(co-1,length(linestyles));
                h(np) = plot(x,y(:,co),'color',colors{np},'linestyle',linestyles{ls});
                hold on;
            end
            for co = 1:size(sdfs.extras)
                if ~isempty(sdfs.extras{co})
                    plot([min(x) max(x)],[sdfs.extras{co}.sdf(sn) sdfs.extras{co}.sdf(sn)],'--','color',colors{np});
                    result.extray(co) = sdfs.extras{co}.sdf(sn);
                    result.extrax(co) = sdfs.extraval(co);
                end
            end
            labels{np} = sprintf('dT = %.0fms',delays(j)/10);
            if size(sdfs.s,2) > 2 & j == bestnv & showpcolor
                GetFigure(labelc);
                hold off;
                pcolor(sdfs.x,sdfs.y,y);
                shading('interp');
                GetFigure(labelb);
            end
            np = np+1;
        end
        if strcmp(Expt.Stimvals.et,'Op') && ~centerRFvals && isfield(Expt.Stimvals,'rOp')
            plot([Expt.Stimvals.rOp Expt.Stimvals.rOp], get(gca,'ylim'),':','color','k');
        end
        if strcmp(Expt.Stimvals.et,'Pp') && ~centerRFvals && isfield(Expt.Stimvals,'rPp')
            plot([Expt.Stimvals.rPp Expt.Stimvals.rPp], get(gca,'ylim'),':','color','k');
        end

        if legendpos < 6
            result.blegend = legend(h,labels,0);
        end
        if getpsych
            title(sprintf('%s Solid = Up',Expt.Header.Name)); %Rd 1, I think
        else
            title(sprintf('%s',Expt.Header.Name));
        end
        if secondorder ==4
            zmax = 0;
            if showpcolor == 2
                subplot(2,2,1);
                hold off;
                zmax(1) = PlotPcolor(sdfs, times, result.delays(cslices(2))-300,interp);
                hold on;
                subplot(2,2,2);
                zmax(2) = PlotPcolor(sdfs, times, result.delays(cslices(1)),interp);
                subplot(2,2,3);
                zmax(3) = PlotPcolor(sdfs, times, result.delays(cslices(2)),interp);
                subplot(2,2,4);
                hold off;
                zmax(4) = PlotPcolor(sdfs, times, result.delays(cslices(3)),interp);
                for j = 1:4;
                    subplot(2,2,j);
                    caxis([0 max(zmax)]);
                    colorbar;
                end
            elseif showpcolor == 1
                hold off;
                PlotPcolor(sdfs, times, result.delays(result.bestdelay),interp);
            else
                subplot(1,1,1);
                h = [];
                labels = {};
                hold off;
                if(autoslice)
                    sn = sampleid(cslices(2));
                else
                    sn = sampleid(slices);
                end
                for gr = 1:length(sn)
                    if length(sn) > 2
                        subplot(2,2,gr);
                    elseif length(sn) > 1
                        subplot(2,1,gr);
                    else
                        subplot(1,1,1);
                    end
                    for pre = 1:size(sdfs.s,3)
                        for k = 1:size(sdfs.s,1)
                            x(k) = sdfs.x(k,1,pre);
                            for co = 1:size(sdfs.s,2)
                                if isempty(sdfs.s{k,co,pre})
                                    y(k,co) = NaN;
                                else
                                    y(k,co) = sdfs.s{k,co,pre}(sn(gr));
                                end
                            end
                        end
                        for co = 1:size(sdfs.s,2)
                            h(pre) = plot(x,y(:,co),'color',colors{pre},'linestyle',linestyles{co});
                            hold on;
                        end
                        labels{pre} = sprintf('pre = %.3f',sdfs.z(1,1,pre));
                    end
                    if(legendpos < 6 & ~isempty(h) && gr == 1)
                        result.legend = legend(h,labels,legendpos);
                    end
                    title(sprintf('At %d',times(sn(gr))));
                end
            end
        end
        

        if profile
            fprintf('Done at %.2f VR %.1f,%.1f,%.1f\n',toc,result.varratio(1:3));
        end
    end
end

if isfield(sdfs,'s')
    result.sdfs = sdfs;
    result.times = times;
    if isfield(Expt.Trials,'Spikes')
    result.delaysamples = sampleid;
    end
end


for j = 1:length(extraidx)
    tx(extraidx{j}) = -(1000+j);
end
for j = 1:length(Expt.Trials)
    Expt.Trials(j).Pd = tx(:,j);
end

result.calctime(2) = now;
if addheader
    result.Header = Expt.Header;
end
% Real end of PlotRevCorAny function


function zmax = PlotPcolor(sdfs, times, t,interp,startt,varargin)

if nargin < 5
    startt = 10;
end
refsdf = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'refsdf',6)
        j = j+1;
        refsdf = varargin{j};
    end
    j = j+1;
end
if isempty(t)   %% Plot first order resp w.r.t t
for j = 1:size(sdfs.x,1)
    for k = startt:length(times)
        ik = k - startt+1;
        X(j,ik) = times(k);
        if sdfs.x(j,1,1) > -1000
            Y(j,ik) = sdfs.x(j,1,1);
        else
            Y(j,ik) = NaN;
        end
        
        if length(refsdf)
            Z(j,ik) = sdfs.s{j,1,1}(k)-refsdf(k);
        else
            Z(j,ik) = sdfs.s{j,1,1}(k);
        end
    end
end
else
id = find(times == t);
for j = 1:size(sdfs.x,1)
    for k = 1:size(sdfs.x,3)
        X(j,k) = sdfs.x(j,1,k);
        Y(j,k) = sdfs.z(j,1,k);
        if id > length(sdfs.s{j,1,k})
            Z(j,k) = NaN;
        elseif length(refsdf)
            Z(j,k) = sdfs.s{j,1,k}(id)-refsdf(id);
        else
            Z(j,k) = sdfs.s{j,1,k}(id);
        end
    end
end
end
[X,Y,Z] = fillpmesh(X,Y,Z);
pcolor(X,Y,Z);
colormap('hot');
if interp
    shading('interp');
else
    shading('flat');
end

title(sprintf('At %d ms',t/10));
zmax = max(caxis);

function [Expt, tidx,ns] = FindTrials(Expt, x, y, tx, ty, good, isextra, varargin)

secondorder =  0;
loop = 0;
extra = 0;
gval = 1;
extray = 0;
j =1;
while j < nargin -6
    if strncmpi(varargin{j},'secondorder',4)
        j = j+1;
        loop = varargin{j}(2);
        secondorder = varargin{j}(1);
        if ismember(secondorder, [4 5 6])
            xval = varargin{j}(3);
        end
    elseif strncmpi(varargin{j},'xvals',4)
        j = j+1;
        xvals = varargin{j};
    elseif strncmpi(varargin{j},'RespDir',5)
        j = j+1;
        gval = varargin{j};
    elseif strncmpi(varargin{j},'extra',4)
        if strncmpi(varargin{j},'extrabyy',8)
            extray = 1;
        end
        j = j+1;
        extra = varargin{j};

    end
    j = j+1;
end

tidx = [];
ns = 0;

for j = 1:length(Expt.Trials)
    if extra
        if extray %% only get extras that match for Y value
            idx = find(isextra(:,j) == extra & good(:,j) == gval & y == ty(:,j));
        else
            idx = find(isextra(:,j) == extra & good(:,j) == gval);
        end
        if secondorder == 1
            if loop == 1
                fidx = find(ty(:,j) == 1);
            else
                fidx = find(ty(:,j) == -1);
            end
            idx = intersect(idx,fidx+1);
        end
    else
% if secondorder ==1, build RC using only frames preceded by
% negative (fy == -1) one where the preceding frame was positive (fy==1)
% i.e. for AC expts, look at resp to Corr when preceded by Corr or AntiCorr
        if secondorder == 1
            if loop == 1
                fidx = find(ty(:,j) == 1);
            else
                fidx = find(ty(:,j) == -1);
            end
            idx = find(tx(:,j) == x & ty(:,j) == y & good(:,j) == gval);
            idx = intersect(idx,fidx+1);
% if secondorder == 2, build 2 RCs, one where the preceding frame had
% x < 0 and y ==1 ie crossed correlated disparity.
        elseif secondorder == 2
            if loop == 1
                fidx = find(tx(:,j)  < 0 & ty(:,j) ==1);
            else
                fidx = find(tx(:,j) > 0 & ty(:,j) ==1);
            end
            idx = find(tx(:,j) == x & ty(:,j) == y & good(:,j) == gval);
            idx = intersect(idx,fidx+1);
% if secondorder == 3, build RC only where the preceding frame had
% x < current (loop = 1) or x > current (loop = 2)
        elseif secondorder == 3
            if loop == 1
                fidx = find(tx(:,j)  < x);
            elseif loop == 2
                fidx = find(tx(:,j) > x);
            elseif loop == 3
                fidx = find(tx(:,j) == 0);
            end
            idx = find(tx(:,j) == x & ty(:,j) == y & good(:,j) == gval);
            idx = intersect(idx,fidx+1);
% if secondorder == 4, build RC for all preceding x vals;
% if aloop is set to n, build only for precding values xvals[aloop:nloop]
        elseif ismember(secondorder,[4 6])
            fidx = find(tx(:,j) == xval);
            idx = find(tx(:,j) == x & ty(:,j) == y & good(:,j) == gval);
            if ~isempty(idx) & ~isempty(fidx)
                idx = intersect(idx,fidx+1);
            else  %if one in empty, there are no intersections.
                idx = [];
            end

% if secondorder == 5, build RC for all preceding x vals, similar
% to secondorder == 4. But now sort by  second value. So that
% if aloop is set to n, and nloop is only 1, show all responses to 
% stimulus n, with a different line for each preceding stimulus.
        elseif secondorder == 5
            fidx = find(tx(:,j) == x);
            idx = find(tx(:,j) == xval & ty(:,j) == y & good(:,j) == gval);
            idx = intersect(idx,fidx+1);
            idx = find(tx(:,j) == x & ty(:,j) == y & good(:,j) == gval);
        else
            idx = find(tx(:,j) == x & ty(:,j) == y & good(:,j) == gval);
        end
% with full sceond order kernels, include extras in the matrix. 
% Easier than having separate matrices for the extras.
        if secondorder < 4
        idx = idx(find(isextra(idx,j) == 0));
        end
    end
    if ~isempty(idx)
        if ~isfield(Expt.Trials,'TrialStart')
            fprintf('%s not build with -trials, Cannot do RC\n');
            return;
        end
        Expt.Trials(j).Trigger = Expt.Trials(j).Start(idx)' - Expt.Trials(j).TrialStart;
        tidx = [tidx j];
        ns = ns+length(idx);
    end
end



function [Expt, tidx,ns] = NewFindTrials(Expt, x, y, tx, ty, good, isextra, varargin)
%
%Bizarely, this is slower, even though the main step is done as on matrix;
secondorder =  0;
loop = 0;
extra = 0;
gval = 1;
j =1;
while j < nargin -6
    if strncmpi(varargin{j},'secondorder',4)
        j = j+1;
        loop = varargin{j}(2);
        secondorder = varargin{j}(1);
        if ismember(secondorder, [4 5])
            xval = varargin{j}(3);
        end
    elseif strncmpi(varargin{j},'xvals',4)
        j = j+1;
        xvals = varargin{j};
    elseif strncmpi(varargin{j},'RespDir',5)
        j = j+1;
        gval = varargin{j};
    elseif strncmpi(varargin{j},'extra',4)
        j = j+1;
        extra = varargin{j};
    end
    j = j+1;
end

tidx = [];
ns = 0;

[fidx, tidx] = find(tx == x & ty == y & good == gval & isextra == 0);
ns = length(tidx);
for j = 1:length(Expt.Trials)
    id = find(tidx == j);
    Expt.Trials(j).Trigger = Expt.Trials(j).Start(fidx(id))' - Expt.Trials(j).TrialStart;
end


