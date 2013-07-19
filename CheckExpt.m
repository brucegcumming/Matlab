function [isi, details] = CheckExpt(Expt, varargin)
%isi = CheckExpt(Expt)
%CheckExpt calculates total # spikes, trial duration, and ISIS for an expt
%isi(1) is mean isi, isi(2) is median.
plotdur = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'plotdur',7)
        plotdur = 1;
    end
    j = j+1;
end
if iscell(Expt) & ischar(Expt{1})
    for j = 1:length(Expt)
        [isi{j}, details{j}] = CheckExpt(Expt{j});
    end
elseif iscell(Expt) & isfield(Expt{1},'Trials') % a cell array of expts
    for j = 1:length(Expt)
        [isi{j}, details{j}] = CheckExpt(Expt{j});
    end
elseif iscell(Expt)
    details = Expt; %plot up existing data
end
if ischar(Expt)
    Expt = LoadExpt(Expt);
end

if iscell(Expt)
    if plotdur
        for j = 1:length(details);
            durs(j,1) = details{j}.tdur + details{j}.ntrials .*2000; %want 100ms each side anyway.
            durs(j,2) = details{j}.TotalDur(1);
        end
    id = find(durs(:,1) > 0 & durs(:,2) > 0);
    plot(durs(id,1),durs(id,2),'o');
    [a,b] = fit_bothsubj2error(durs(id,1),durs(id,2));
    refline(b);
    title(sprintf('Slope %.1f',b));
    end
    return;
end
delay = 500;
nerr = 0;
rptframes = [];
for j = 1:length(Expt.Trials)
    ends(j) = Expt.Trials(j).End(end);
    starts(j) = Expt.Trials(j).Start(1);    
    durs(j) = Expt.Trials(j).End(end) - Expt.Trials(j).Start(1);
    spks(j) = length(find(Expt.Trials(j).Spikes > delay & Expt.Trials(j).Spikes < durs(j)+delay));
    if isfield(Expt.Trials,'rptframes')
        rptframes(j) = length(Expt.Trials(j).rptframes);
    end
end
isis = starts(2:end) - ends(1:end-1);
isi(1) = mode(isis);
isi(2) = median(isis);

if ~isempty(rptframes)
    id = find(rptframes > 0);
    details.nrpt = length(id);
    details.meanrpt = mean(rptframes(id));
end

details.tdur = sum(durs);
details.TotalDur(1) = Expt.Trials(end).End(end) - Expt.Trials(1).Start(1);
details.ntrials = length(Expt.Trials);
if isfield(Expt.Header,'BlockStart') && length(Expt.Trials) > 1
    if Expt.Trials(end).Trial > Expt.Header.BlockStart(end)+1
        Expt.Header.BlockStart = [Expt.Header.BlockStart Expt.Trials(end).Trial];
    end
    for j = 1:length(Expt.Header.BlockStart)-1
        ta = find([Expt.Trials.Trial] == Expt.Header.BlockStart(j));
        tb = find([Expt.Trials.Trial] == Expt.Header.BlockStart(j+1))-1;
        edur(j) = Expt.Trials(tb).End(end) - Expt.Trials(ta).Start(1);
    end
    details.TotalDur(1) = sum(edur);
details.TotalDur(2) = Expt.Trials(end).End(end) - Expt.Trials(1).Start(1);
else
details.TotalDur = Expt.Trials(end).End(end) - Expt.Trials(1).Start(1);  
end
       wrongsacs = [];
       wrongids = [];
if isfield(Expt.Trials,'RespDir') & isfield(Expt.Trials,'Saccades')

    for j = 1:length(Expt.Trials)
        len = Expt.Trials(j).End - Expt.Trials(j).Start;
        idx = find([Expt.Trials(j).Saccades.start] > len & abs([Expt.Trials(j).Saccades.size]) > 1);
        if isempty(idx) 
                if Expt.Trials(j).RespDir ~= 0
                rdirs(j) = NaN;
                else
                   rdirs(j) = 0;
                end
        else
               
        if Expt.Trials(j).RespDir * sin(Expt.Trials(j).Saccades(idx(1)).dir) < 0
            wrongsacs = [wrongsacs j];
            wrongids = [wrongids idx(1)];
        else
            Saccade(j) = Expt.Trials(j).Saccades(idx(1));
        end
        if(Expt.Trials(j).RespDir == 0)
                rdirs(j) = 0;
        else
rdirs(j) =  Expt.Trials(j).Saccades(idx(1)).dir;
        end
        end
    end
elseif isfield(Expt.Trials,'RespDir')
    nerr = nerr+1;
    details.errs{nerr} = 'Respdir but no saccaedes';
    fprintf('%s %s\n',Expt.Header.Name,details.errs{nerr});
end
fprintf('%d spikes mean rate %.1f\n',sum(spks),10000 * mean(spks./durs));
if ~isempty(wrongsacs)
    
for j = 1:length(wrongsacs)
    trial = Expt.Trials(wrongsacs(j));
    sac = trial.Saccades(wrongids(j));
    fprintf('Trial %.0f Resp %.0f Sac %.1f at %.0f\n',wrongsacs(j),trial.RespDir,sac.dir,sac.start-(trial.End-trial.Start));
end
end