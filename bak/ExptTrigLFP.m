function [avg, details] = ExptTrigLFP(Expt,LFP, varargin)
% avg = ExptTrigLFP(Expt,LFP, varargin)
% build a spike triggerd LFP for an Expt file

blocks = [Expt.Header.BlockStart Expt.Trials(end).Trial];
tid = [];
nsplit = 0;
argon = {};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'split',4)
        j = j+1;
        splittype = varargin{j};
        j = j+1;
        splitvals = varargin{j};
        nsplit = length(splitvals);
    else
        argon = {argon{:} varargin{j}};
    end
    j = j+1;
end

if isfield(Expt,'probes') %% copmbine across probes, 
else
for j = 1:length(Expt.Header.Clusters)
    if size(Expt.Header.Clusters{j},2) < Expt.Header.probe | isempty(Expt.Header.Clusters{j}{1,Expt.Header.probe})
        autocut(j) = 1;
    elseif ~isfield(Expt.Header.Clusters{j}{1,Expt.Header.probe},'autocut')
        autocut(j) = 0;
    else
        autocut(j) = Expt.Header.Clusters{j}{1,Expt.Header.probe}.autocut;
    end
    if autocut(j) == 0
        tid = [tid blocks(j):blocks(j+1)];
    end
end
end

LFP = CheckLFP(LFP,'fix');

if isempty(tid) %% all autocut. use all
    tid = 1:length(Expt.Trials);
    lid = find(ismember([LFP.Trials.Trial],[Expt.Trials.Trial]));
else
    etrials = tid;
    tid = find(ismember([Expt.Trials.Trial],etrials));
    lid = find(ismember([LFP.Trials.Trial],[Expt.Trials(tid).Trial]));
end

if isempty(lid)
    fprintf('No LFP Trials for %s\n',Expt.Header.Name);
    avg = [];
    details.tid = tid;
    details.lfpid = [LFP.Trials.Trial];
    return;
else
    LFP.Trials = LFP.Trials(lid);
end
for j = 1:length(LFP.Trials)
    LFP.Trials(j).Spikes = Expt.Trials(tid(j)).Spikes;
end
if isfield(LFP.Trials,'RespDir') & isfield(LFP.Trials,'Dc')
    tid = find([LFP.Trials.Dc] < 0.01 & [LFP.Trials.RespDir] < 0);
    [avg(:,:,1), details(1)] = SpTrigLFP(LFP.Trials(tid),20000,1./(LFP.Header.LFPsamplerate.*10000),100,argon{:});
    tid = find([LFP.Trials.Dc] < 0.01 & [LFP.Trials.RespDir] > 0);
    [avg(:,:,2), details(2)] = SpTrigLFP(LFP.Trials(tid),20000,1./(LFP.Header.LFPsamplerate.*10000),100,argon{:});
elseif isfield(LFP.Trials,'RespDir') & isfield(LFP.Trials,'ob')
    tid = find([LFP.Trials.ob] > 120 & [LFP.Trials.RespDir] < 0);
    if length(tid)
    [avg(:,:,1), details(1)] = SpTrigLFP(LFP.Trials(tid),20000,1./(LFP.Header.LFPsamplerate.*10000),100,argon{:});
    else
        details(1).npk = 0;
        avg = [];
    end
    tid = find([LFP.Trials.ob] > 120 & [LFP.Trials.RespDir] > 0);
    if length(tid)
    [avg(:,:,2), details(2)] = SpTrigLFP(LFP.Trials(tid),20000,1./(LFP.Header.LFPsamplerate.*10000),100,argon{:});
    else
        details(2).nspk = 0;
    end
elseif nsplit > 0 & isfield(LFP.Trials,splittype);
    nc = 1;
    for k = 1:length(splitvals)
        if isnan(splitvals(k)) && isfield(LFP.Trials,'st')
    tid = find(mean([LFP.Trials.st],1) < 1);            
        else
    tid = find(mean([LFP.Trials.(splittype)],1) == splitvals(k));
        end
    if length(tid)
    [avg(:,:,nc), details(nc)] = SpTrigLFP(LFP.Trials(tid),20000,1./(LFP.Header.LFPsamplerate.*10000),100,argon{:});
    nc = nc+1;
    end
    end
else
[avg, details] = SpTrigLFP(LFP.Trials,20000,1./(LFP.Header.LFPsamplerate.*10000),100,argon{:});
end

        