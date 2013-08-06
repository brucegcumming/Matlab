function [pp, details, Expt] = ExptPsych(Expt, varargin)
%function [pp, details] = ExptPsych(Expt, varargin)
%  Fit pyschometric funtion to an Expt file that has responses.
%
%ExptPsych(Expt, 'trials', list)
%   uses only Expt.trials(list)
%ExptPsych(Expt, 'trialid', list)
%   uses only Expt.Trials.Trial in list
% if Expt is a cell array, each element is plotted separately.
%N.B. to plot a subgroup, be sure the argument is still a cell array, not a
%set of cells
%ExptPsych(Expts(2:3),....)   is what works
%
%ExptPsych(Expts(2:3),'sum')  
%                      combines these expts into one.
%ExptPsych(Expt,'smooth',[a b])
%                       plots response rates for each stim over time,
%                       smoothed with boxcar width a. And % correct at max
%                       signal, smoothed with b.
%ExptPsych(Expt,'pcolor',[n b])
%                       plots functions for binned groups of n trails
%                       as a pseudocolor plot. If b is given, the bin is
%                       moved this distane between samples
%ExptPsych(Expt,'pcolorfit',[n b])
%                       fits psf to each group.
%
%Except in recognized formats, the expt1 variable is used as the metameter.
%If expt2 is set, the usual default is to collapse across these. To force
%separate plotting, use ..,'collapse',[0 0 0]

            
labelb = 'Sequence';

fitpredict = 0;

if iscell(Expt)
    if length(varargin) > 0 && strncmpi(varargin{1},'sum',3)
        E = CombineExpts(Expt);
        [pp, details] = ExptPsych(E,varargin{2:end});

       return;
    else
    labels = {};
    h = [];
    k = 1;
    for j = 1:length(Expt)
        if sum([Expt{j}.Trials.RespDir] ~= 0) > 20

        [pp{k}, details{k}] = ExptPsych(Expt{j},'color',j-1,varargin{:});
        if sum(isfield(details{k},{'labels' 'handles'})) ==2 && length(details{k}.labels) == length(details{k}.handles)
           labels = {labels{:} details{k}.labels{:}};
           h = [h details{k}.handles];
        end
        drawnow;
        k = k+1;
        end
    end
    if length(labels) > 1
        mylegend(h,labels);
    end
    return;
    end
end
type{1} = 'Dc';
type{2} = 'ori';
type{3} = '';
showplot = 1;
psfargs = {};
pp = [];
details = [];
smoothw = [0 0];
plotreward = 0;
plotsequence = 0;
plotfit = 0;
shown = 0;
pntrials = 0;
ptag = 'Psych Pcolor';
tag = '';
maintag = [];
nofit = 0;
plottype = 1;
minreps = 4;
signals = [];
collapse = [0 1 0]; %collapse exp2 by default. Becauce?? 
coloroff = 0;
verbose = 0;
mintrials = 10;
predictchoice = 0;
predictcrit = NaN;
prednoise = 0;
labelth = 0;
labelbias = 0;
show.th = 0;
show.bias = 0;
choiceonly = 1;
legendlabels = {};

typestr = '';
explabel = ''; % fot the title
if strcmp(Expt.Stimvals.et,'or') && strcmp(Expt.Stimvals.e2,'ob')
    type{2} = 'or';
    for j = 1:length(Expt.Trials)
        if Expt.Trials(j).ob < 0
            Expt.Trials(j).or = Expt.Trials(j).or - 90;
        end
        Expt.Trials(j).cv = sd2cv(abs(Expt.Trials(j).ob));
    end
    type{1} = 'cv';
    typestr = 'ORBW';
elseif strcmp(Expt.Stimvals.et,'dx')
    if strcmp(Expt.Stimvals.e2,'Dc')
        type{2} = 'dx';
        type{1} = Expt.Stimvals.e2;
        typestr = 'DCDP';
        collapse = [0 0 0];
    else
    if sum(strcmp(Expt.Stimvals.e2,{'rb' 'or' 'mixac' 'cz' 'co' 'Dc'}))
        collapse = [0 0 0];
    end
    type{1} = 'dx';
    type{2} = Expt.Stimvals.e2;
    typestr = 'DT';
    if sum(strcmp(Expt.Stimvals.e2,{'rb' 'or' 'mixac' 'cz' 'co'}))
        collapse = [0 0 0];
    end
    end
elseif strcmp(Expt.Stimvals.et,'TwoCylDisp')
    type{1} = 'dx';
    type{2} = 'rd';
    for j = 1:length(Expt.Trials)
        Expt.Trials(j).rd = Expt.Trials(j).dx-Expt.Trials(j).bd;
    end
    typestr = 'DT';
    collapse = [0 0 1];
    if Expt.Stimvals.Bc < 1
    explabel = sprintf('Or%.0f, Back co %.2f',Expt.Stimvals.or,Expt.Stimvals.Bc);
    else
    explabel = sprintf('Or%.0f, Back or%.0f',Expt.Stimvals.or,Expt.Stimvals.bo);
    end
else
    type{1} = Expt.Stimvals.et;
    type{2} = Expt.Stimvals.e2;
    collapse = [0 0 0];
end


if isfield(Expt.Trials,'rwsum') %online data
    explabel = [explabel sprintf('rw%.1f',Expt.Trials(end).rwsum)];
end
if strcmp(type{1},'Dc') && strcmp(type{2},'or')
    type{2} = 'ori';
end

if strmatch(type{2},'ori')  & ~isfield(Expt.Trials,'ori')
    for j = 1:length(Expt.Trials)
        Expt.Trials(j).ori = Expt.Trials(j).or(end);
    end
end

fillet = { [] [] []};
if isfield(Expt.Stimvals,'et') && ~isfield(Expt.Trials,Expt.Stimvals.et)
    fillet{1} = Expt.Stimvals.et;
end
if isfield(Expt.Stimvals,'e2') && ~isfield(Expt.Trials,Expt.Stimvals.e2)
    fillet{2} = Expt.Stimvals.e2;
end
if isfield(Expt.Stimvals,'e3') && ~isfield(Expt.Trials,Expt.Stimvals.e3)
    fillet{3} = Expt.Stimvals.e3;
end

for j = 1:length(Expt.Trials)
        if isfield(Expt.Trials,'Start')
        starts(j) = Expt.Trials(j).Start(1);
        else
        starts(j) = j;
        end
        if ~isfield(Expt.Trials,'Trial') || isempty(Expt.Trials(j).Trial)
            Expt.Trials(j).Trial = j;
        end
        Expt.Trials(j).e0 = 0;
        if isfield(Expt.Trials,'exvals')
            for k = 1:length(fillet)
                if ~isempty(fillet{k})
                Expt.Trials(j).(fillet{k}) = Expt.Trials(j).exvals(k);
                end
            end
        end
end
    
    
if length(varargin) & strcmp(varargin{end},'byrw')
    [nrw, rws] = Counts([Expt.Trials.rw]);
    rws = rws(rws > 0);
    if length(rws) < 5
        hold off;
        for j = 1:length(rws)
            [pp, psf(j)] = ExptPsych(Expt,varargin{1:end-1},'rwsz',rws(j),'color',j);
            labels{j} = sprintf('rw%.2f',rws(j));
            hold on;
        end
        mylegend([psf.handles],labels);
    else
    hold off; [pp, a] = ExptPsych(Expt,varargin{1:end-1},'rwmin',0.14);
    hold on;  [pp, b] = ExptPsych(Expt,varargin{1:end-1},'rwmax',0.081,'color',2);
    if size(a.fit.fit) > 2
    c.shifts = [a.fit.fit(3)-a.fit.fit(1) b.fit.fit(3)-b.fit.fit(1)];
    c.slopes = [a.fit.fit(2) b.fit.fit(2)];
    details = {a b c};
    t = get(get(gca,'title'),'string');
    title(sprintf('%s/%.4f Shift %.4f/%.4f',t,c.slopes(1),c.shifts(2),c.shifts(1)));
    end
    end
    return;
end


%Expt 3 types to be plotted separately by default
if isfield(Expt.Stimvals,'e3') && sum(strcmp(Expt.Stimvals.e3,{'tr' 'nf' 'uS'}))
    type{3} = Expt.Stimvals.e3;
end

if isfield(Expt.Trials,'tr') %some guided responses included
    if strcmp(type{2},'e0')
        type{2} = 'tr';
        collapse(2) = 0;
    end
end

j = 1;
while j <= length(varargin)
   if strncmpi(varargin{j},'nmin',3)
       j = j+1;
       minreps = varargin{j};
       psfargs = {psfargs{:} 'nmin' varargin{j}};
   elseif strncmpi(varargin{j},'color',5)
       j = j+1;
       coloroff = varargin{j};
   elseif strncmpi(varargin{j},'collapse',5)
       j = j+1;
       collapse = varargin{j};
   elseif strncmpi(varargin{j},'fitpredict',7)
       fitpredict = 1;
       predictchoice = 1;
   elseif strncmpi(varargin{j},'labelth',7)
       show.th = 1;
   elseif strncmpi(varargin{j},'labelbias',7)
       show.bias = 1;
   elseif strncmpi(varargin{j},'legendlabels',7)
       j = j+1;
       legendlabels = varargin{j};
   elseif strncmpi(varargin{j},'mintrials',5)
       j = j+1;
       mintrials = varargin{j};
   elseif strncmpi(varargin{j},'nofit',3)
       nofit = 1;
   elseif strncmpi(varargin{j},'noplot',3)
       showplot = 0;
   elseif strncmpi(varargin{j},'pcolor',3)
       if strncmpi(varargin{j},'pcolorfit',8)
           nofit = 0;
           plotfit = 1;
       else
           nofit = 1;
       end
       j = j+1;
       pntrials = varargin{j};
       if length(pntrials) == 1
           pntrials(2) = round(pntrials(1)/10);
       end
   elseif strncmpi(varargin{j},'predicted',6)
       predictchoice = 1;
       if length(varargin) > j && isnumeric(varargin{j+1})
           j = j+1;
           predictcrit = varargin{j};
       end
   elseif strncmpi(varargin{j},'prednoise',6)
       predictchoice = 1;
       j = j+1;
       prednoise = varargin{j};
   elseif strncmpi(varargin{j},'reward',3)
       plotreward = 1;
   elseif strncmpi(varargin{j},'sequence',3)
       plotsequence = 1;
   elseif strncmpi(varargin{j},'shown',5)
       psfargs = {psfargs{:} 'shown'};
       shown = 1;
   elseif strncmpi(varargin{j},'show',3)
       j = j+1;
       show.(varargin{j}) = 1;
   elseif strncmpi(varargin{j},'signals',3)
       j = j+1;
       signals = varargin{j};
   elseif strncmpi(varargin{j},'smooth',3)
       j = j+1;
       smoothw =  varargin{j};
       if length(smoothw) == 1
           smoothw(2) = smoothw(1);
       end
   elseif strncmpi(varargin{j},'rwsz',4)
       j = j+1;
       id = find([Expt.Trials.rw] == varargin{j});
       Expt.Trials = Expt.Trials(id);
   elseif strncmpi(varargin{j},'rwmax',5)
       j = j+1;
       id = find([Expt.Trials.rw] <= varargin{j});
       Expt.Trials = Expt.Trials(id);
   elseif strncmpi(varargin{j},'rwmin',5)
       j = j+1;
       id = find([Expt.Trials.rw] >= varargin{j});
       Expt.Trials = Expt.Trials(id);
   elseif strncmpi(varargin{j},'tag',3)
       j = j+1;
       tag = varargin{j};
   elseif strncmpi(varargin{j},'type2',5)
       j = j+1;
       type{2} = varargin{j};
   elseif strncmpi(varargin{j},'type3',5)
       j = j+1;
       type{3} = varargin{j};
   elseif strncmpi(varargin{j},'trialids',8)
       j = j+1;
       usetrials = varargin{j};
       id = find(ismember([Expt.Trials.Trials],usetrials));
       Expt.Trials = Expt.Trials(id);
   elseif strncmpi(varargin{j},'trials',3)
       j = j+1;
       Expt.Trials = Expt.Trials(varargin{j});
   elseif strncmpi(varargin{j},'verbose',3)
       verbose = 1;
   elseif strncmpi(varargin{j},'xmax',3)
       j = j+1;
       psfargs = {psfargs{:} 'xmax' varargin{j}};
   elseif strncmpi(varargin{j},'xmin',3)
       j = j+1;
       psfargs = {psfargs{:} 'xmin' varargin{j}};
   end
    j = j+1;
end

if ~isfield(Expt.Header,'expname')
    if isfield(Expt.Header,'name')
        Expt.Header.expname = Expt.Header.name;
    else
        Expt.Header.expname = Expt2Name(Expt);
    end
end

if ~isfield(Expt.Trials,type{1}) || (~isfield(Expt.Trials, type{2}) && length(type{2}))
    fprintf('%s Doesn''t have values for Psych\n',Expt.Header.expname);
    return;
end


nt = sum(abs([Expt.Trials.RespDir]) > 0);
if nt < mintrials
    if isfield(Expt.Header,'expname')
        fprintf('%s too few response trials\n',Expt.Header.expname);
    elseif isfield(Expt.Header,'Name')
        fprintf('%s too few response trials\n',Expt.Header.Name);
    end
    return;
end

if choiceonly
id = find(abs([Expt.Trials.RespDir]) ==1); %only use choice trials by 
Expt.Trials = Expt.Trials(id);
end

sigvar = type{1};

predchoices = [];
if predictchoice && isfield(Expt.Trials,'pkresp')
    crits = [];
    id = find(abs([Expt.Trials.RespDir]) ==1);
    if ~isnan(predictcrit)
        crit = predictcrit;
        if length(predictcrit) > 1
            crits = predictcrit;
        end
    else
        zid = find([Expt.Trials.ob] == 130 & abs([Expt.Trials.RespDir]) == 1);
        crit = prctile([Expt.Trials(zid).pkresp],50);
        predictcrit = crit; %%save to use as guess for fitting
        prednoise = std([Expt.Trials(zid).pkresp])./3;
    end
    if isempty(crits)
        crits = ones(size(id)).* crit;
    end
    noise = randn(size(id)) .* prednoise;
    pkresps = [Expt.Trials(id).pkresp];
    for j = 1:length(id)
        if Expt.Trials(id(j)).pkresp +noise(j)> crits(j)
            predchoices(id(j)) = 1;
  %          Expt.Trials(id(j)).RespDir = 1;
        elseif Expt.Trials(id(j)).pkresp +noise(j) < crits(j)
%            Expt.Trials(id(j)).RespDir = -1;
            predchoices(id(j)) = -1;
        end
    end
end

if ~isempty(tag)
    GetFigure(tag);
end

if strcmp(typestr,'ORBW') & isfield(Expt.Trials,'se')
    id = find(abs([Expt.Trials.cv]) < 0.01 & ismember([Expt.Trials.RespDir],[-1 1]));
    [a,b,c] = CalcConsistency(Expt.Trials(id));
    details.consistency = a;
    details.seedrpts = c;
    details.seeds = b;
    if ~isfield(Expt.Header,'BlockStart')
        Expt.Header.BlockStart = 1;
    end
    BlockEnd = [Expt.Header.BlockStart(2:end)-1  Expt.Trials(end).Trial];
    
    for j = 1:length(Expt.Header.BlockStart)
        id = find(abs([Expt.Trials.cv]) < 0.01 & ...
            ismember([Expt.Trials.RespDir],[-1 1])...
            & [Expt.Trials.Trial] >= Expt.Header.BlockStart(j) ...
            & [Expt.Trials.Trial] <= BlockEnd(j));
        [a,b,c] = CalcConsistency(Expt.Trials(id));
        details.consistencies(j) = a;
        a = sum([Expt.Trials(id).RespDir] == 1)./length(id);
        details.biases(j) = a;
        details.n(j) = length(id);
        details.nrpt(j) = sum(c(c > 1));
        id = find([Expt.Trials.Trial] >= Expt.Header.BlockStart(j) ...
            & [Expt.Trials.Trial] <= BlockEnd(j));
        [Expt.Trials(id).blockid] =deal(j);
    end
end

if strcmp(Expt.Stimvals.e3,'Us')
    if isempty(type{2})
        type{2} = 'e0';
        vals{2} = zeros(size([Expt.Trials.(type{1})]));
    end
        type{3} = 'uStim';
else
    vals{3} = zeros(size([Expt.Trials.(type{1})]));
    uvals{3} = [0];
end

if length(type) < 3 && isfield(Expt.Trials,'tr') && length(unique([Expt.Trials.tr])) > 1
    type{3} = 'tr';
    uvals{3} = unique([Expt.Trials.tr]);
end

for j = 1:length(type)
    if isfield(Expt.Trials,type{j})
    vals{j} = [Expt.Trials.(type{j})];
    uvals{j} = unique(vals{j});
    end
end

if nargout > 2
        for j = 1:length(Expt.Trials)
            Expt.Trials(j).signal = abs(Expt.Trials(j).(sigvar)) .* Expt.Trials(j).rwdir;
        end
end

if verbose
    fprintf('%d Expts...',length(type));
end
[counts, range] = hist(abs(vals{1}),unique(abs(vals{1})));
if plotreward
    scores = [Expt.Trials.RespDir].* [Expt.Trials.rwdir];
    plot(scores,'o');
    hold on;
    plot([Expt.Trials.rw]./max([Expt.Trials.rw]*1.1),'ro-');
    return
elseif plotsequence
    respdir = [Expt.Trials.RespDir];
    signal = abs([Expt.Trials.(sigvar)]) .* [Expt.Trials.rwdir];
    for j = 1:length(Expt.Trials)
        if length(uvals{2}) > 1
            sgn = sign(Expt.Trials(j).(type{2}) - mean(uvals{2}));
        else
            sgn = 1;
        end
        Expt.Trials(j).signal = signal(j);
    end
    signal = [Expt.Trials.signal];
    signals = unique([Expt.Trials.signal]);
    scores = [Expt.Trials.RespDir].* [Expt.Trials.rwdir];
    rtype = 2.5 +[scores + respdir/2];
    for j = 1:length(signals)
    %wrong L, wrong R, correct L, correct R, 
    rid = find(signal == signals(j) & respdir > 0);
    lid = find(signal == signals(j) & respdir < 0);
    rid = rid(rid > 1);
    lid = lid(lid > 1);
    if isempty(rid) || isempty(lid)
        loseshift(j) = NaN;
        winstay(j) = NaN;
    else
    nt(j) = length(rid)+length(lid);
    k(:,1) = histc(rtype(rid-1),[1:4])./length(rid);
    k(:,2) = histc(rtype(lid-1),[1:4])./length(lid);
    nwrong(j) = sum(sum(k(1:2,:)));
    ncorrect(j) = sum(sum(k(3:4,:)));
    loseshift(j) = (k(1,1)+k(2,2)-k(1,2)-k(2,1))./nwrong(j);
    winstay(j) = (k(3,1)+k(4,2)-k(3,2)-k(4,1))./ncorrect(j);
    end
    end
    [a,zid] = min(abs(signals));
    details.ntrials = nt;
    details.signals = signals;
    details.loseshift = cat(1,loseshift, nwrong);
    details.winstay = cat(1,winstay, ncorrect);
    F = gcf;
    GetFigure(labelb);
    hold off; 
    plot(details.signals,loseshift,'o-');
    hold on;
    plot(details.signals,winstay,'ro-');
    legend('Lose-Shift','Win-Stay');
    title(sprintf('%d Trials. LoseShift %.1f/%.2f, Winstay %.1f/%.2f',nt(zid),loseshift(zid).*100,nwrong(zid)/2,winstay(zid).*100,ncorrect(zid)/2));
    figure(F);
elseif pntrials(1) > 0
    trials = [Expt.Trials.Trial];
    k = 1;
    step = pntrials(2);
    args = {};
    if nofit
        args = {args{:} 'nofit'}
    end
    pp= {};
    details = {};
    for j = 1:step:(length(trials)-pntrials(1))
        [pp{k}, details{k}] = ExptPsych(Expt,'trials',[j:j+pntrials(1)-1],args{:});
        lens(k) = length(pp{k});
        tn(k) = j; 
        k = k+1;
    end
    nex = length(lens);
    if plotfit
        for j = 1:nex
            X(j,:) = details{j}.fit.fitx;
            Y(j,:) = tn(j) .* ones(size(X(j,:)));
            Z(j,:) = details{j}.fit.fity;
            slopes(j) = details{j}.fit.fit(1);
            bias(j) = details{j}.fit.fit(2);
        end
    else
    X = ones(k+1,max(lens)+1) .* NaN;
    Y = ones(k+1,max(lens)+1) .* NaN;
    Z = ones(k+1,max(lens)+1) .* NaN;
    for j = 1:nex
        [xs, xid] = sort([pp{j}.x]);
        for k = 1:lens(j)
            X(j,k) = pp{j}(xid(k)).x;
            Y(j,k) = tn(j);
            Z(j,k) = pp{j}(xid(k)).resp./pp{j}(xid(k)).n;
        end
            X(j,k+1) = pp{j}(xid(k)).x + mean(diff([pp{j}.x]));
            Y(j,k+1) = tn(j);
            Z(j,k+1) = pp{j}(xid(k)).resp./pp{j}(xid(k)).n;
    end
    end
    f = gcf;
    GetFigure(ptag);
    hold off; 
    pcolor(X',Y',Z');
    shading('flat');
    colorbar;
    caxis([0 1]);
    colormap('hot');
    figure(f);
    return;
elseif smoothw
   hold off;
   colors = mycolors;
    if strcmp(type{1},'cv')
        signs = vals{2} - mean(uvals{2});
        signvals{1} = vals{1} .* sign(signs);
    else 
        signvals{1} = vals{1};
    end
   for j = 1:length(uvals{1})
       for k = 1:length(uvals{2})
           id = find(vals{1} == uvals{1}(j) & vals{2} == uvals{2}(k) & [Expt.Trials.RespDir] ~= 0);
           if length(id) > smoothw(1)
            [x,y] = xysmooth(starts(id)', [Expt.Trials(id).RespDir]',smoothw(1));
            y = (y + 1)./2;
           h(j) = plot(x,y,'o-','color',colors{j});
           labels{j} = sprintf('%.3f',uvals{1}(j));
           hold on;
           end
       end
   end
   id = find(ismember([Expt.Trials.Trial],Expt.Header.BlockStart));
   if length(id) > 1
   for j = 1:length(id) 
       j = Expt.Trials(id(j)).Start(1);
       plot([j j],[0 1],'k:');
   end
   end
   id = find(range > max(range)/2);
   pid = find(signvals{1} >= min(range(id)) & abs([Expt.Trials.RespDir]) > 0);
   nid = find(signvals{1} <= -min(range(id)) & abs([Expt.Trials.RespDir]) > 0);
   apid = find(signvals{1} > 0  & abs([Expt.Trials.RespDir]) > 0);
   anid = find(signvals{1} < 0 & abs([Expt.Trials.RespDir]) > 0);
   if mean([Expt.Trials(nid).RespDir]) < 0 
       scores(apid) = ([Expt.Trials(apid).RespDir] +1)/2;
       scores(anid) = (1-[Expt.Trials(anid).RespDir])/2;
   elseif mean([Expt.Trials(nid).RespDir]) > 0 
       scores(apid) = (1- [Expt.Trials(apid).RespDir])/2;
       scores(anid) = ([Expt.Trials(anid).RespDir]+1)/2;
   end
   for j = 1:length(uvals{1})
       id = find(vals{1} == uvals{1}(j) & [Expt.Trials.RespDir] ~= 0);
      [x,y] = xysmooth(starts(id)', scores(id),smoothw(2));
       plot(x,y,'color',colors{j},'linewidth',2);
   end

   id = union(pid,nid);
   [x,y] = xysmooth(starts(id)', scores(id),smoothw(2));
   j = length(h)+1;
   h(j) = plot(x,y,'k','linewidth',2);
   labels{j} = sprintf('>%s',num2str(max(range)/2));
   mylegend(h,labels);
   if isfield(Expt.Trials,'rb')
       plot([Expt.Trials.Start],(0.9-[Expt.Trials.rb])./1.8);
   end
   id = find(abs([Expt.Trials.RespDir]) > 0);
   [x,y] = xysmooth(starts(id)', scores(id),smoothw(2));
   h(j) = plot(x,y,'k:','linewidth',2);
   return;
end

signs = sign(uvals{2} - mean(uvals{2}));
labels = {};
np = 1;
for j = 1:length(uvals{1})
    for k = 1:length(uvals{2});
        for m = 1:length(uvals{3});
            if length(uvals{3} > 1)
                 si = find(vals{1} == uvals{1}(j) & vals{2} == uvals{2}(k) & vals{3} == uvals{3}(m));
            else
                si = find(vals{1} == uvals{1}(j) & vals{2} == uvals{2}(k));
            end
        rup = sum([Expt.Trials(si).RespDir] < 0);
        rdown = sum([Expt.Trials(si).RespDir] > 0);
        pp(np).n = rup+rdown;
        pp(np).resp = rup;
        pp(np).y = uvals{2}(k);
        pp(np).z = uvals{3}(m);
        pp(np).tid = si;
        tids{np} = si;
        if strcmp(type{2},'e0')
            pp(np).x = uvals{1}(j);
            pp(np).expno = m;
        elseif strcmp(type{2},'Id') && strcmp(type{1},'bd')
            pp(np).x = uvals{1}(j);
            if uvals{2}(k) > 1
                pp(np).expno = 4;
            elseif uvals{2}(k) > 0
                pp(np).expno = 1;
            elseif uvals{2}(k) < 0
                pp(np).expno = 2;
            else
                pp(np).expno = 3;
            end
        elseif strcmp(type{2},'Id') 
            pp(np).x = uvals{1}(j);
            if abs(pp(np).x) > 0  || uvals{2}(k) == 0 %% should test for mD
                pp(np).expno = 1;
                labels{1} = 'ID0';
            elseif uvals{2}(k) > 1
                pp(np).expno = 4;
                pp(np+1) = pp(np);
                np = np+1;
                pp(np).expno = 1;
                labels{4} = 'IDx';                
            elseif uvals{2}(k) > 0
                pp(np).expno = 2+(m)*length(uvals{3});
                pp(np+1) = pp(np);
                np = np+1;
                pp(np).expno = 1;
                labels{2} = 'ID+';
            else
                pp(np).expno = 3+(m)*length(uvals{3});
                pp(np+1) = pp(np);
                np = np+1;
                pp(np).expno = 1;
                labels{3} = 'ID-';
            end
        elseif sum(strcmp(type{2},{'ori' 'dx'}))
            pp(np).x = uvals{1}(j) .* signs(k);
            pp(np).expno = m;            
        else
% don't need signs for most expts. Only if expt 1 is Dc/ob
           if strcmp(type{1},'cv')
            pp(np).x = uvals{1}(j) .* signs(k);

           else
            pp(np).x = uvals{1}(j);
           end
            if collapse(2)
                pp(np).expno = m;
            else
                pp(np).expno = k;
                if strcmp(typestr,'ORBW') & uvals{2}(k) >= 30
                    pp(np).expno = find(uvals{2} == uvals{2}(k)-90);
                end
            end
        end
        if length(uvals{3}) > 1
            pp(np).expno = m;
        end
        np = np+1;
        end
    end
end


if type{1} == 'cv'
    id = find(abs([pp.x]) < 0.01);
    if length(id) > 1
        pp(id(1)).n = sum([pp(id).n]);
        pp(id(1)).resp = sum([pp(id).resp]);
        pp(id(1)).tid = cat(2,pp(id).tid);
        [pp(id(2:end)).n] = deal(0);
        [pp(id(2:end)).resp] = deal(0);
        [pp(id(2:end)).tid] = deal([]);
    end
    ors = unique(vals{2});
    for j = 1:length(ors)/2
        for m = 1:length(uvals{3})
        id = find([pp.y] == ors(j) & [pp.z] == uvals{3}(m));
        ei = j + (m-1) * length(uvals{3});
        [pp(id).expno] = deal(ei);
        id = find([pp.y] == ors(j)+90 & [pp.z] == uvals{3}(m));
        [pp(id).expno] = deal(ei);
        end
    end
end
if verbose
    fprintf('%d conditions',length(pp));
end
for j = unique([pp.expno])
    for k = 1:length(uvals{1})
    id = find([pp.x] == uvals{1}(k) & [pp.expno] == j);
    if length(id) > 1
    pp(id(1)).n = sum([pp(id).n]);
    pp(id(1)).resp = sum([pp(id).resp]);
    uid = setdiff(1:length(pp),id(2:end));
    pp = pp(uid);
    end
    end
end

if nofit
    return;
end

if length(predchoices)
    for j = 1:length(pp)
        ppp(j) = pp(j);
        ppp(j).resp = sum(predchoices(tids{j}) == -1);
    end
end
if length(uvals{3}) > 1 && sum(strcmp(type{3},{'tr' 'nf'})) == 0
fit = fitpsf(pp,'twomeans');
if showplot
    ph = fitpsf(fit.data,'twomeans', fit,psfargs{:},'coff',coloroff,'showfit');
    if verbose
        fprintf('Fit ');
    end
    h = ph.fith; %handles to fit lines
    if strmatch(type{3},'Ustim')
    labels{1} = 'NoStim';
    labels{2} = 'uStim';
    else
        labels{1} = sprintf('%s=%.2f',type{3},uvals{3}(1));
        labels{2} = sprintf('%s=%.2f',type{3},uvals{3}(2));
    end
    set(gca,'Ylim',[0 1]);
    title(sprintf('Bias %.4f,%.4f (=%.4fp < %.3f), Slope %.4f',fit.fit(1),fit.fit(3),fit.fit(1)-fit.fit(3),fit.pval,abs(fit.fit(2))));
    xlabel(type{1});
    text(min([pp.x]),-0.1,0,sprintf('%s %.1f',type{2},uvals{2}(1)));
    text(max([pp.x]),-0.1,0,sprintf('%s %.1f',type{2},uvals{2}(end)),...
        'HorizontalAlignment','right');
end
else
    yv = unique([pp.expno]);

    cols = mycolors;
    for j = 1:length(yv)
        icolor = j+coloroff;
        id = find([pp.expno] == yv(j) & [pp.n] >= minreps);
        fits{j} = fitpsf(pp,'expno',yv(j),psfargs{:});
        if showplot & isfield(fits{j},'data')
            id= find([pp.expno] == yv(j));
            fit = fits{j};
            plot([fit.data.x],[fit.data.p],'o')
            Expt.Stimvals.th = abs(fit.fit(2));
            Expt.Stimvals.bias = fit.fit(1);
%            plot([pp(id).x],[pp(id).resp]./[pp(id).n],'o');
            if length(id) > 1
                ph = fitpsf(fit.data,'showfit',fit,psfargs{:},'color',cols{j+coloroff});
                fit.fith = ph.fith;
                h(j) = ph.fith;
                if j > length(labels)
                    if length(legendlabels) >= j
                        labels{j} = legendlabels{j};
                    elseif length(uvals{3}) > 1
                        labels{j} = [sprintf('%s=',type{3}) sprintf('%.3f ',unique([pp(id).z]))];
                    else
                        labels{j} = [sprintf('%s=',type{2}) sprintf('%.3f ',unique([pp(id).y]))];
                    end
                    f = fields(show);
                    for k = 1:length(f)
                        if show.(f{k}) && isfield(Expt.Stimvals,f{k})
                            labels{j} = [labels{j} f{k} '=' num2str(Expt.Stimvals.(f{k}))];
                        end
                    end
                    
                end
            else
                h(j) = plot([pp(id).x],[pp(id).resp]./[pp(id).n],'s','color',cols{j+coloroff});
                if shown
                    text([pp(id).x],[pp(id).resp]./[pp(id).n],num2str(pp(id).n));
                end
            end
            n = [pp.n];
            pp = pp(n > 0);
            set(gca,'Ylim',[0 1]);
            if shown
                explabel = [num2str(sum(n)) 'Trials ' explabel];
            end
            title(sprintf('Bias %.1f, Slope %.4f %s',fit.fit(1),abs(fit.fit(2)),explabel));
            xlabel(type{1});
            if strmatch(type{2},{'ob'})
                text(min([pp.x]),-0.1,0,sprintf('%s %.1f',type{2},uvals{2}(1)));
                text(max([pp.x]),-0.1,0,sprintf('%s %.1f',type{2},uvals{2}(end)),...
                    'HorizontalAlignment','right');
            end
        if length(predchoices)
            if fitpredict
                blkid = [Expt.Trials.blockid];
                if isfield(details,'consistency')
                    [pp.consistency] = deal(details.consistency);
                end
                options = optimset('MaxFunEvals',100000,'maxiter',1000,'display','off');
                guess(1) = mean(pkresps);
          %      guess(1)=0;
                guess(2) = std(pkresps)./20;
              %  guess(2) = 0;
              guess(1) = predictcrit;
              guess(2) = prednoise;
              guess(3) = guess(1);
                setappdata(gcf,'params',[]);
                [ssd, ppp] = MinimizeCrit(guess, pp, pkresps,blkid);
                [fittedparams,fval,exitflag, output] = fminsearch(@MinimizeCrit,guess,options,pp, pkresps,blkid);
                [ssd, ppp] = MinimizeCrit(fittedparams, pp, pkresps,blkid);
                details.fitpp = ppp;
                details.predfit.params = fittedparams;
                params  = getappdata(gcf,'params');
                
            end

            if showplot
                pfits{j} = fitpsf(ppp,'showfit','expno',yv(j),psfargs{:},'color',cols{j+coloroff+1});
            else
                pfits{j} = fitpsf(ppp,'showfit','expno',yv(j),psfargs{:});
            end
            details.predfit.fits = pfits;
        end
        elseif length(id) == 1
            if length(uvals{3}) > 1
                icolor = find([pp(id).z] == uvals{3}) + coloroff;
            else
                icolor = j+coloroff;
            end
            plot(pp(id).x,pp(id).resp./pp(id).n,'s','color',cols{icolor});
            hold on;
            l = text(pp(id).x,pp(id).resp./pp(id).n,sprintf('%s=%.2f n= %d',type{2},pp(id).y,pp(id).n));
        end

    end
end
    if verbose
    fprintf('%d Expts',length(yv));
    end
if exist('fits','var')
details.fits = fits;
end
if exist('fit','var')
details.fit = fit;
details.handles = h;
end
id = find([pp.n] > 0);
minn = median([pp(id).n]) * 0.8;
id = find([pp.n] > minn);
perf = [pp(id).resp]./[pp(id).n];
details.maxperf = mean([1-min(perf) max(perf)]);
if ~exist('ph','var')
    return;
end
if showplot && length(id)
    details.fit.fitx = ph.fitx;
    details.fit.fity = ph.fity;
end
details.labels = labels;
if isfield(Expt.Trials,'rw')
[a,b] = Counts([Expt.Trials.rw]);
details.rws = b;
details.rwn = a;
end
if strcmp(type{3},'nf') && length(h) == length(uvals{3}) && length(h) > 1
    labels = {};
    for j = 1:length(h)
        labels{j} = sprintf('%s=%.1f',type{3},uvals{3}(j)); 
    end
end

if length(labels) == length(h);
    id = find(h > 0 & ishandle(h));
    if details.fit.fit(2) < 0
        legend(h(id),labels{id});
    else
        legend(h(id),labels{id},'Location','SouthEast');
    end
end

if strcmp(Expt.Stimvals.e3,'Us') && Expt.Header.rc
    GetFigure('PK');
    hold off;
    if isempty(signals)
    xv = sort(unique([Expt.Trials.(type{1})]));
    signals = xv(xv < 0.1);
    end
    id = find([Expt.Trials.uStim] > 0);
    [K(1,:), a] = ExptPRC(Expt,signals, id);
    id = find([Expt.Trials.uStim] == 0);
    rK(1,:) = (a.ratios' * sum(a.nstims,2))/sum(a.nstims(:));
    [K(2,:), b] = ExptPRC(Expt, signals, id);
    rK(2,:) = (b.ratios' * sum(b.nstims,2))/sum(b.nstims(:));
    id= find(b.rcvals >=0);
    plot(b.rcvals(id),K(:,id)');
    hold on;
    plot(b.rcvals(id),diff(K(:,id)),'r');
    if length(signals) > 1
    scale = max(abs(K(:)))./(max(abs(rK(:)-0.5)));
    plot(b.rcvals(id),(rK(:,id)-0.5) .* scale,':');
    end
    legend('Stim','non','diff');
elseif Expt.Header.rc
    [K,b] = ExptPRC(Expt, 0, []);
    GetFigure('PK');
    hold off;
    id= find(b.rcvals >=0);
    plot(b.rcvals(id),K(id))
end    


function [consistency, ses, srpts] = CalcConsistency(Trials)
    ses = unique([Trials.se]);
    samedir = 0;
    diffdir = 0;
    srpts = [];
    for j = 1:length(ses)
        sid = find([Trials.se] == ses(j));
        srpts(j) = length(sid);
        for k = 2:length(sid) 
            if Trials(sid(k)).RespDir == Trials(sid(k-1)).RespDir
                samedir = samedir+1;
            else
                diffdir = diffdir+1;
            end
        end
    end
    consistency = samedir./(samedir+diffdir);
    

function [kernel, details]  = ExptPRC(Expt, signals, id, varargin)
kernel = [];
details = [];

if isempty(id)
    id = 1:length(Expt.Trials);
end
for j = 1:length(id)
    lens(j) = length(Expt.Trials(id(j)).or);
end
len = prctile(lens,90);
sid = find(lens < len);
for j = sid;
    Expt.Trials(id(j)).or(lens(j)+1:len) = NaN;
end

R.dcs = [Expt.Trials(id).Dc];
signs = sign([Expt.Trials(id).ori] - mean([Expt.Trials(id).ori]));
R.dcs = R.dcs .*signs;
signals = union(signals,signals * -1); %need to do + and - signals separaely
R.resp = [Expt.Trials(id).RespDir];
R.seq = [Expt.Trials(id).or]';
id = find(R.seq > -999);
R.seq(id) = mod(R.seq(id),180);
[kernel, details] = CalcPRC(R,signals, varargin{:});
id = find(sum(details.nstims,2) > 0);
ns = details.nstims(id,:);
kernel = sum(kernel(id,:),1)./sum(ns(:));

function [kernel, details] = CalcPRC(R, signals, varargin)

colors = {'b' 'r' 'm'};
plotmeans = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'colors',5)
        j = j+1;
        colors = varargin{j};
    end
    j = j+1;
end

ival = 0;
for j = 1:length(signals)
    signal = signals(j);
upid = find(R.resp == 1 & R.dcs == signal);
dnid = find(R.resp == -1 & R.dcs == signal);
aid = find(R.dcs == signal);
sv = unique(R.seq(:));
sv = sv(find(~isnan(sv)));
svid = 1:length(sv);
dnseq = R.seq(dnid,:);
upseq = R.seq(upid,:);
uk = hist(upseq(:),sv(svid))./length(upid);
dk = hist(dnseq(:),sv(svid))./length(dnid);

if (length(dnid) == 0 || length(upid) == 0)
    kernel(j,:)  = zeros(size(sv));
    nstim(j,:) = [0 0];
    details.ratios(j,:) = zeros(size(sv));
else
vcounts = sum(R.seq(aid,:)' == ival);
ak = (uk.*length(upid) + dk .* length(dnid))./(length(upid) + length(dnid));
pid = find(uk+dk > 0);
nstim(j,:) = [length(upid) length(dnid)];
kernel(j,:) = (dk-uk).*sum(nstim(j,:));
details.ratios(j,:) = dk./(uk+dk);
dks(j,:) = dk;
uks(j,:) = uk;
end
end  % end for signals


%what is this?  finds all R.seq == 0, then counts them up
%? a check that p(blank) was working?  
vvals = unique(vcounts);
if sum(vvals) > 0  %not emppty
for j = 1:length(vvals)
    id = find(vcounts == vvals(j));
    presp(j) = mean(R.resp(aid(id)));
end
end
mk = sum(kernel,1)./sum(nstim(:));
alln = sum(dks,1)+sum(uks,1);

id = find(alln > max(alln)/(3 + 12 * abs(signal)));
if plotmeans == 0
    id = id(find(~ismember(sv(svid(id)),[-1 1 2]))); % remove extras
end
details.dks = dks;
details.uks = uks;
details.rcvals = sv;
details.nstims = nstim;
id = id(find(~ismember(sv(svid(id)),[-1 1 2]))); % remove extras
if length(id) > 1
[a, details.sineamp] = famp(sv(svid(id)), mk(id), 1/180);
else 
    details.sineamp = NaN;
end

function p = ChoiceP(crit, signal, noisesd)


g = 0.5 * (1 + erf((signal-crit)./sqrt(2)./noisesd));
p = g;
id = find(g < 0.5);
g(id) = 1-g(id);

function [SSD, ppp]  = MinimizeCrit(x,pp,pkresps, blkid)

crit = ones(size(pkresps)) .* x(1);
simulate = 0;
if simulate
    nr=100;
    if length(x) == 3
        blks = unique(blkid);
        for b = blks
            id = find(blkid ==b);
            crit(id(floor(length(id/2)):end)) = x(3);
        end
    end
    crit = repmat(crit,nr,1);
    pkresps = repmat(pkresps,nr,1);
    pkresps = pkresps+ randn(size(pkresps)).*x(2);
else
    nr=1;
    pkp = 1 - ChoiceP(crit,pkresps,x(2));
    id = find(pkp < 0.5);
    psame = pkp;
    psame(id) = 1-psame(id);
end

id = find(pkresps < crit);
predchoices = ones(size(pkresps));
predchoices(id) = -1;
for j = 1:length(pp)
    pp(j).p = pp(j).resp./pp(j).n;
    if abs(pp(j).x) < 0.001 && pp(j).n> 0
        zstim = j;
    end
end
for j = 1:length(pp)
    ppp(j) = pp(j);
    ppp(j).resp = sum(sum(predchoices(:,pp(j).tid(:)) == -1))./nr;
    if simulate
        ppp(j).p = ppp(j).resp./ppp(j).n;
    else
        ppp(j).p = mean(pkp(pp(j).tid(:)));
        ppp(j).resp = ppp(j).p .* ppp(j).n;
        ppp(j).consistency = mean(psame(pp(j).tid(:)));
    end
    
end

if simulate
[ppp.consistency] = deal(0);

if nr > 1
    for k = 1:length(ppp)
    tid = ppp(k).tid;
    for j = 1:length(tid)
        same(j) = sum(predchoices(:,tid(j)) ==1)./nr;
        if same(j) < 0.5;
            same(j) = 1 - same(j);
        end
    end
    ppp(k).consistency = mean(same);
    end
end
end
SSD = sum(([ppp.resp]-[pp.resp]).^2);
if length(x) == 3
    SSD = SSD + ((ppp(zstim).consistency-pp(zstim).consistency) .* pp(zstim).n).^2;
end
params = getappdata(gcf,'params');
params = cat(1, params, [x SSD]);
setappdata(gcf,'params',params);


