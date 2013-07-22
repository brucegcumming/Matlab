function [result] = PlotRateSequence(Expts, varargin)

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
   
   j = 1;
   while j <= length(varargin) 
       if strncmpi(varargin{j},'bytime',5)
           xstyle = 'time';
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
       end
       j = j+1;
   end
if isstruct(Expts) && isfield(Expts,'Trials');
    Expt = Expts;
    clear Expts;
    Expts{1} = Expt;
end
   
   for j = 1:length(Expts)
    Expt = Expts{j};
    dur = mean([Expt.Trials.dur])./10000;
    rates = [Expt.Trials.count]./dur;
    counts = [counts rates];
    result.meanrates(j) = mean(rates);
    tn = [tn Expt.Trials.Trial];
    ends(j) = tn(end);
            E = Expts{j};
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
            AllBlocks(j) =  E.Trials(1).TrialStart;
        else
        AllTrials = cat(2,AllTrials,[E.Trials.Trial]);
        AllBlocks(j) =  E.Trials(1).Trial;
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
        dy = 0;
        scale = 1;
    end
    result.meanrates = (result.meanrates .* scale) + dy;
    result.rates = (scale.*counts)+dy;
    if isempty(callback)
        h  = plot(AllTrials, (scale.*counts)+dy,'o','color',color);
    else
        h  = plot(AllTrials, (scale.*counts)+dy,'o','color',color,'buttondownfcn',callback);
    end
    result.AllBlocks = AllBlocks;
    result.handle = h;
    
   
    function HitPoint(a,b,t)
        
        fprintf('Trial %d\n',t);
        
    