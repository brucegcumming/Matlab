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
Expt = LoadGridLFP(Expt, varargon{:});
Expt = LoadClusterInfo(Expt);
Expt = FixLFPTrials(Expt,varargon{:});
if fixspike
    Expt = FixLFPSpike(Expt, varargon{:});
end


