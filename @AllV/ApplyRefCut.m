function ApplyRefCut(a,b, varargin)

DATA = GetDataFromFig(a);
DATA.autorefine = 0;
uselastcluster = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{1},'lastcluster',6)
        uselastcluster = 1;
    elseif strncmpi(varargin{1},'refine')
        DATA.autorefine = 1;
    end
    j = j+1;
end


p = DATA.probelist(DATA.probe(1));

if uselastcluster
    Clusters = AllV.mygetappdata(DATA,'LastClusters');
    DATA.cluster = Clusters{p};
    fprintf('Using Expt %d Clusters for probe %d\n',p,Clusters{p}.exptno);
    DATA.userefcluster = 0;

else
    if ~isfield(DATA,'RefClusters') || length(DATA.RefClusters) < p
        DATA.RefClusters = LoadRefClusters(DATA);
    end
    DATA.cluster = DATA.RefClusters{p};
    DATA.userefcluster = 1;
    fprintf('Using RefClusters for probe %d\n',p);
end

DATA.cluster.manual = 4;
DATA.cluster.strictscaling = 1; %re-adjust template scaling so SDs are the eame
DATA = AllV.ClassifyAll(DATA, 1);