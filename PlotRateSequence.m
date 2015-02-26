function [result] = PlotRateSequence(Expts, varargin)
%PlotRateSequence(Expts, ...) plots the rates for each trial in a set of expts
%(different blocks for one cell)
%PlotRateSequence(AllExpt) Plots sequences for each cell in AllExpt

tn = [];
counts = [];
AllTrials = [];
AllIds = [];
result = [];
AllBlocks = [];
xstyle = 'trial';
normalize = 0;
   offset = 0;
   color = 'b';
   callback = {};
plotcells = 0;
checkcounts = 0;
state.tag = 'MeanRates';
smoothw = 0;

   j = 1;
   while j <= length(varargin) 
       if strncmpi(varargin{j},'bytime',5)
           xstyle = 'time';
       elseif strncmpi(varargin{j},'cells',5) %Expt expt{} is a different cell
           plotcells = 1;
       elseif strncmpi(varargin{j},'check',5)
           checkcounts = 1;
       elseif strncmpi(varargin{j},'callback',8)
           j = j+1; 
           callback = varargin{j};
       elseif strncmpi(varargin{j},'color',5)
           j = j+1; 
           color = varargin{j};
       elseif strncmpi(varargin{j},'normalize',5)
           normalize = 1;
       elseif strncmpi(varargin{j},'offset',5)
           j = j+1;
           offset = varargin{j};
       elseif strncmpi(varargin{j},'smooth',5)
           j = j+1;
           smoothw = varargin{j};
       end
       j = j+1;
   end
   
%?? this causes infinite recusion with iscell() below
%why was it here
if isstruct(Expts) && isfield(Expts,'Trials') && 0
    Expt = Expts;
    clear Expts;
    Expts{1} = Expt;
end
if isstruct(Expts) && isfield(Expts,'allrates')
    GetFigure(state.tag);
    result = CheckCounts(Expts, state, varargin{:});
    return;
end
if iscellstr(Expts)
    names = Expts;
    for j = 1:length(names)
        Expt = LoadExpt(names{j});
        if isempty(Expt)
            result{j}.name = names{j};
            AddError(result{j},'Empty Expt in %s',names{j});
        else
            result{j} = PlotRateSequence(Expt,varargin{:});
        end
    end
    return;
elseif iscell(Expts) && ~isempty(CellToMat(Expts,'maxtrial')) %array of results
    for j = 1:length(Expts)
        result{j} = PlotRateSequence(Expts{j},varargin{:});
    end
    return;
elseif ischar(Expts)
    
    name = Expts;
    Expt = LoadExpt(name);
    result = PlotRateSequence(Expt,varargin{:});
    result.Expt = Expt;
return;
end


if isstruct(Expts) && isfield(Expts,'Spikes');
    A = Expts;
    clear Expts;
    Expts = All2Expt(A,'all');
    filename = A.Expt.Header.loadname;
    plotcells = 1;
    result.maxtrial = max([A.Expt.Trials.Trial]);
end

cellnumber = [];
if plotcells
    colors = mycolors;
    y = 0;
    hold off;
    for j = 1:length(Expts)
        if isfield(Expts{j},'plotres')
            E = Expts{j}.plotres.Data;
        else
            E = Expts{j};
        end
        if ~isempty(E.Trials)
        if normalize
            x{j} = PlotRateSequence(E,'color',colors{j}, 'offset', (j-1),'normalize');
            line(get(gca,'xlim'),[j-1 j-1],'color',colors{j});
        else
            x{j} = PlotRateSequence(E,'color',colors{j}, 'offset', y);
            line(get(gca,'xlim'),[y y],'color',colors{j});
            y = y+max([E.Trials.count]);
        end
        if isfield(E,'trialsused')
            x{j}.uset = E.trialsused;
        else
            x{j}.uset = 1:length(x{j}.rates);
        end
        x{j}.ids = [E.Trials.id];
        allids{j} = x{j}.times(:);
        cellnumber(j) = E.Header.cellnumber;
        hold on;
        end
    end
    if checkcounts
        ntrials = length(unique(cat(1,allids{:})));
        alltrials = CellToMat(x,'times'); %trial #, usually
        triallist = unique(alltrials);
        triallist = triallist(triallist > 0);
        allrates = ones(length(x),ntrials).*NaN;
        for j = 1:length(x)
            allrates(j,x{j}.uset) = x{j}.rates - x{j}.offset;
            ids(x{j}.uset) = x{j}.ids;
        end
        result.allrates = allrates;

        a = sum(allrates ==0);
        id = find(a > length(x)/4);
        if ~isempty(id)
            result = AddError(result,'%s: Zero rates on %d-%d probes Trials %s\n',GetEval(E,'shortname'),min(a(id)),max(a),sprintf(' %d',ids(id)));
        end
        result.trials = triallist;
        result.ids = ids;
        result.missingid = ids(id);
        result.missing = id;
        result.name = GetName(E,'path');
        result.checktime = now;
        triallist = triallist(triallist > 0);
        result.cellnumber = cellnumber;
        for j = 1:length(triallist)
            id = find(alltrials == triallist(j));
            meanrate(j) = mean(allrates(id));
            meanvar(j) = var(allrates(id));
        end
        GetFigure(state.tag);
        hold off;
        plot(triallist,meanrate,'k-');
        hold on;
        plot(triallist,meanvar,'r-');
    end
    return;
end

%This plots a list of Expts for one cell.
   for j = 1:length(Expts)
       if iscell(Expts)
           Expt = Expts{j};
       else
           Expt = Expts(j);
       end
    dur = mean([Expt.Trials.dur])./10000;
    rates = [Expt.Trials.count]./dur;
    counts = [counts rates];
    result.meanrates(j) = mean(rates);
    tn = [tn Expt.Trials.Trial];
    ends(j) = tn(end);
            E = Expt;
        AllIds = cat(2,AllIds,[E.Trials.id]);
        if strcmp(xstyle,'time')
            if isfield(E.Header,'timeoffset')
                toff = E.Header.timeoffset;
                if isfield(E.Header,'timeadjust')
                    toff = E.Header.timeoffset-E.Header.timeadjust;
                end
            else
                toff = 0;
            end
            AllTrials = cat(2,AllTrials,toff+([E.Trials.TrialStart]./10000));
            if isempty(E.Trials)
                AllBlocks(j) =  E.Header.trange(1);
            else
                AllBlocks(j) =  E.Trials(1).TrialStart;
            end
        else
        AllTrials = cat(2,AllTrials,[E.Trials.Trial]);
        AllBlocks(j) =  E.Trials(1).Trial;
        if isfield(E.Header,'cellnumber')
            cellnumber(j) = E.Header.cellnumber;
        end
        end
   end
   result.times = AllTrials;
    [AllTrials, id] = unique(AllTrials);
    AllIds = AllIds(id);
    counts = counts(id);
    dy = 0;
    if normalize == 1
        scale = 1./mean(counts);
        dy = offset;
    elseif normalize == 2
        dy= mean(counts) .* offset;
        scale = 1;
    else
        dy = offset;
        scale = 1;
    end
    result.meanrates = (result.meanrates .* scale) + dy;
    result.rates = (scale.*counts)+dy;
    result.offset = dy;
    result.scale = scale;
    result.cellnumbers = cellnumber;
    result.meancount = mean(counts);
    if isempty(callback)
        h  = plot(AllTrials, (scale.*counts)+dy,'o','color',color);
    else
        h  = plot(AllTrials, (scale.*counts)+dy,'o','color',color,'buttondownfcn',callback);
    end
    if smoothw
        hold on;
        h  = plot(AllTrials, smooth((scale.*counts)+dy,smoothw),'-','color',color);
    end
    result.AllBlocks = AllBlocks;
    result.handle = h;
    
function x = CheckCounts(x, state, varargin)
    
    cellsonly = 0;
    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'cellsonly',4)
            cellsonly = 1;
        end
        j = j+1;
    end
    showsd = 0;
    a = sum(x.allrates ==0);
    id = find(a > size(x.allrates,1)/4);
    if ~isempty(id)
        fprintf('%s: Zero rates on %d-%d probes Columns %s\n',GetName(x),min(a(id)),max(a),sprintf(' %d',x.ids(id)));
    end
    %x.missingid = ids(id);
    x.missing = id;
    x.checktime = now;
    if cellsonly && isfield(x,'cellnumber')
        cid = find(x.cellnumber > 0);
    else
        cid = 1:size(x.allrates,1);
    end
    triallist = x.trials;
    for j = 1:size(x.allrates,2)
        meanrate(j) = nanmean(x.allrates(:,j));
        meanvar(j) = nanvar(x.allrates(:,j));
    end
    GetFigure(state.tag);
    hold off;
    id = find(~isnan(meanrate));
    plot(triallist(id),meanrate(id),'k-');
    if showsd
        hold on;
        plot(triallist(id),meanvar(id),'r-');
    end
    x.zrange = minmax(zscore(meanrate(id)));

    
    function HitPoint(a,b,t)
        
        fprintf('Trial %d\n',t);
        
    