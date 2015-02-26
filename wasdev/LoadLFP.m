function Expt = LoadLFP(Expt, varargin)
% Expt = LoadLFP(Expt, ...)
%Load LFP data into Trials 
%Pads data with NaN so that length of LFP is the same in each trial, and this
%matches with Expt.Header.preperiod (if defined).
% Expt = LoadLFP(Expt, 'zeropad') pads with zeros instead of NaN

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

if iscell(Expt)
    for j = 1:length(Expt)
        if isfield(Expt{j},'Header') && isfield(Expt{j},'Trials')
        Expt{j} = LoadLFP(Expt{j},varargin{:});
        end
    end
    return;
end

if strcmp(Expt.Header.DataType,'Spike2')
    Expt= LoadSpike2LFP(Expt,'reload','fixshort',10);
    if align
        Expt = FixLFPTrials(Expt,varargon{:});
    end
else
    Expt = LoadGridLFP(Expt, varargon{:});
    Expt = LoadClusterInfo(Expt);
    Expt = FixLFPTrials(Expt,varargon{:});
if fixspike
    Expt = FixLFPSpike(Expt, varargon{:});
end
end


