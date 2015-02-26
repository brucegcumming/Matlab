function [Expt, Summary] = LoadLFP(Expt, varargin)
% Expt = LoadLFP(Expt, ...)
%Load LFP data into Trials 
%Pads data with NaN so that length of LFP is the same in each trial, and this
%matches with Expt.Header.preperiod (if defined).
% Expt = LoadLFP(Expt, 'zeropad') pads with zeros instead of NaN
% Expt = LoadLFP(Expt, 'double') makes Expt.Trials.LFP a double (default is int16). 
% [Expt, Summary] = LoadLFP(Expt) returns a structure with informatoin
% about errors in the LFP files (missing data etc
% see also FullV2LFP, FixLFPTrials


fixspike = 0;
align = 1;
varargon = {};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'dummy',5)
        for k = 1:length(Expt.Trials)
            Expt.Trials(k).LFP = rand(1100,24);
            Expt.Trials(k).ftime = Expt.Trials(k).Start(1);
        end
        Expt.Header.CRsamplerate = 0.1;
        Expt.Header.LFPtimes = [1:1100]./1000;
        return;
    elseif strncmpi(varargin{j},'fixspike',6)
        fixspike = 1;
    elseif strncmpi(varargin{j},'noalign',6)
        align = 0;
    else
        varargon = {varargon{:} varargin{j}};
    end
    j = j+1;
end


if iscellstr(Expt)
    names = Expt;
    for j = 1:length(names)
        if isdir(names{j})
            E = ReadExptDir(names{j});
            [Expts{j}, Summary{j}] = LoadLFP(E, varargin{:});
        end
    end
    return;
elseif iscell(Expt)
    for j = 1:length(Expt)
        if isfield(Expt{j},'Header') && isfield(Expt{j},'Trials')
            [Expt{j}, Summary{j}] = LoadLFP(Expt{j},varargin{:});
        end
    end
    for j = 1:length(Expt)
        if isfield(Expt{j},'Header') && isfield(Expt{j}.Header,'LFPmissing') && (Expt{j}.Header.LFPmissing > 0 || Expt{j}.Header.LFPclipped > 0)
            fprintf('Expt%d %d/%d Trials missing LFP %d/%d Start clipped\n',GetExptNumber(Expt{j}),...
                Expt{j}.Header.LFPmissing,length(Expt{j}.Trials),...
                Expt{j}.Header.LFPclipped,length(Expt{j}.Trials));
        end
    end
    return;
end

Summary.gotlfp = 0;
if isfield(Expt,'Spikes') && length(Expt.Header) > 1 %an allexpt 
    [Expt.Expt, Summary] = LoadLFP(Expt.Expt);
    Summary.Name = GetName(Expt);
elseif strcmp(Expt.Header.DataType,'Spike2')
    [Expt, details] = LoadSpike2LFP(Expt,'reload','fixshort',10,varargon{:});
    if align
        Expt = FixLFPTrials(Expt,varargon{:});
    end
    Summary.Name = GetName(Expt);
    Summary.exptno = GetExptNumber(Expt);
    Summary.ntrials = length(Expt.Trials);
    Summary = CopyFields(Summary,Expt.Header,'LFPMissing', 'LFPclipped');
    Summary = CopyFields(Summary,Expt,'LFPerrs','LFPerrdata');    
    Summary = CopyFields(Summary,details,'gotlfp');
else
    [Expt, loaded] = LoadGridLFP(Expt, varargon{:});
    Summary.Name = GetName(Expt);
    Summary.gotlfp = loaded;
    Expt = LoadClusterInfo(Expt);
    if sum(loaded) == 0
        Expt = AddError(Expt,'No LFP Files available. Run FullV2LFP');
    else
        if sum(loaded ==0)
            Expt = AddError(Expt,'%d LFP Files missing.  Run FullV2LFP',sum(loaded==0));
        end
        Expt = FixLFPTrials(Expt,varargon{:});
        if fixspike
            Expt = FixLFPSpike(Expt, varargon{:});
        end
    end
end


