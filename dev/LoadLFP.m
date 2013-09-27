function Expt = LoadLFP(Expt, varargin)
% Expt = LoadLFP(Expt, ...)
%Load LFP data into Trials 

fixspike = 0;
varargon = {};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'fixspike',6)
        fixspike = 1;
    else
        varargon = {varargon{:} varargin{j}};
    end
    j = j+1;
end

if iscell(Expt)
    for j = 1:length(Expt)
        Expt{j} = LoadLFP(Expt{j},varargin{:});
    end
    return;
end

if strcmp(Expt.Header.DataType,'Spike2')
    Expt= LoadSpike2LFP(Expt,'reload','fixshort',10);
else
    Expt = LoadGridLFP(Expt, varargon{:});
    Expt = LoadClusterInfo(Expt);
    Expt = FixLFPTrials(Expt,varargon{:});
if fixspike
    Expt = FixLFPSpike(Expt, varargon{:});
end
end


