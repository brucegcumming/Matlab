function DATA = ReadFinishedDir(DATA, name, varargin)

DataTypes = {};
Trialids = [];
[Expts, Idx] = ReadExptDir(name, varargin{:});
for nexp = 1:length(Expts)
    Expts{nexp}.gui.classified = 0;
    Expts{nexp}.gui.counted = 0;
    Expts{nexp}.gui.clustertype = 0;
    newt = [Expts{nexp}.Trials.Trial];
    Expts{nexp}.gui.firsttrial = newt(1);
    if isfield(Idx{nexp},'errs')
        Expts{nexp}.errs = Idx{nexp}.errs;
    elseif isempty(Idx{nexp})
        cprintf('red','Empty Idx for Expt %d, may have alignement errors\n',nexp);
    end
    Expts{nexp}.gui.ntrials = length(newt);
    if isfield(Idx{nexp},'DataType')
        DataTypes{nexp} = Idx{nexp}.DataType;
    end
    Trialids = [Trialids newt];
    if isfield(Expts{nexp}.Header,'nprobes')
        nprobes(nexp) = Expts{nexp}.Header.nprobes;
    elseif isfield(DATA.ArrayConfig,'id')
        nprobes(nexp) = length(DATA.ArrayConfig.id);
    else
        nprobes(nexp) = 1;
    end
    eids(nexp) = GetExptNumber(Expts{nexp});
end
%
%Empty Expts are NOT included. Suffix list maps included expts to exptno

DATA.suffixlist = eids;
DATA.Expts = Expts;


if ~isempty(DataTypes)
    DATA.DataType = DataTypes{1};
end

DATA.probelist = 1:length(DATA.ArrayConfig.id);
DATA.probenames = cellstr(int2str(DATA.probelist(:)));
for j = DATA.probelist
    DATA.probes(j).probe = j;
end
DATA.probesource = cmb.FindProbeSources(DATA);
DATA = cmb.AddMultiProbeGUI(DATA);
if sum(DATA.probelist < 100) > 1 % count # of recording chans
    DATA.state.includeprobename = 1;
end
%DATA.state.online = 0;
DATA.name = name;
if ~isfield(DATA,'prefix')
    DATA.prefix = name;
end
DATA = cmb.cListExpts(DATA,Expts);
cmb.SetProbeList(DATA);
DATA.exabsid = 1:length(DATA.Expts);
DATA = cmb.CheckCellExptList(DATA);
DATA.AllData.Spikes = [];
DATA.AllData.Trialids = Trialids;
DATA.state.nospikes = 2;
DATA.state.somespikes = 1;
ShowExptErrs(DATA.Expts);

