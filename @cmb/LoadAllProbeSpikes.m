function DATA = LoadAllProbeSpikes(a,b, varargin)
if isstruct(a)
    DATA = a;
else
    DATA = GetDataFromFig(a);
end
useselected = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'select',5)
        useselected = 1;
    end
    j = j+1;
end

if DATA.state.online
    %    DATA = cmb.LoadAllProbesSpikes(DATA);
    %Load spike times for selected experiements ? or do all expts
    %prob best to do all - seelcting the current expts does most of the work
    %anyway????
    args = {};
    if useselected
        args = {args{:} 'expts' DATA.exabsid};
    end
    DATA = cmb.BuildAllTimes(DATA,args{:});
    set(DATA.toplevel,'UserData',DATA);
    return;
end

expid = DATA.currentexpt(1);
times(1) = DATA.Expts{expid}.Trials(1).Start(1)-10000;
times(2) = DATA.Expts{expid}.Trials(end).End(end) + 10000;
DATA.allexp = expid;
DATA.AllData.Spikes = [];
probelist = DATA.probelist;
j = 1;
while j <= length(varargin)
    %        'select' load just selected probes for selected expt range
    if strncmpi(varargin{j},'select',5)
        cmb.AddProbeList(DATA);
        if ~isfield(DATA.plot,'useprobe')
            return;
        end
        probelist = find(DATA.plot.useprobe);
        if isempty(probelist)
            return;
        end
        if sum(probelist == DATA.probe) == 0
            DATA.probe = probelist(1);
        end
        %            DATA.probe = max(probelist);
        eid = DATA.exabsid;
        times(1) = DATA.Expts{eid(1)}.Trials(1).Start(1) - 10000;
        times(2) = DATA.Expts{eid(end)}.Trials(end).End(end) - 10000;
    elseif  strncmpi(varargin{j},'allselect',6) %setlect all
        DATA.plot.useprobe = ones(size(DATA.probelist));
        cmb.AddProbeList(DATA);
        probelist = find(DATA.plot.useprobe);
        eid = DATA.exabsid;
        times(1) = DATA.Expts{eid(1)}.Trials(1).Start(1) - 10000;
        times(2) = DATA.Expts{eid(end)}.Trials(end).End(end) - 10000;
    end
    j = j+1;
end

if isfield(DATA,'AllClusters')
    DATA = cmb.UpdateAllClusters(DATA);
    set(DATA.toplevel,'UserData',DATA);
    return;
end
tic;
DATA = cmb.GetAllProbeFig(DATA);
if ~isfield(DATA,'ptsize')
    DATA.ptsize = 1;
end
if isfield(DATA,'AllClusters') %% Can't have both
    DATA = rmfield(DATA,'AllClusters');
end
for j = 1:length(probelist)
    set(DATA.toplevel,'Name',sprintf('Loading probe %d',probelist(j)));
    drawnow;
    if DATA.state.online == 0
        DATA.AllSpikes{probelist(j)} = cmb.GetProbeFiles(DATA, probelist(j),DATA.subprobe,'trange',times/10000,'nodv');
        %            DATA.AllSpikes{probelist(j)} = cmb.GetProbeFiles(DATA, probelist(j),'trange',times/10000);
        DATA = cmb.LoadClusters(DATA, cmb.ClusterFile(DATA,'probe', probelist(j)),DATA.subprobe,'probe', probelist(j));
    else
        filename = ['C:' DATA.Expts{DATA.currentexpt(1)}.Header.Name];
        if DATA.probelist(j) > 16
            filename = strrep(filename,'/Expt','A/Expt');
        end
        DATA.AllSpikes{probelist(j)} = cmb.GetProbeSpikes(DATA.AllData, filename , DATA.probevars{probelist(j)},[DATA.probe DATA.subprobe]);
        DATA.AllSpikes{probelist(j)}.firstspki = 1;
    end
end
toc
if DATA.state.recut
    tic
    fprintf('Classifying Spikes');
    set(DATA.toplevel,'Name',sprintf('Classifying Spikes'));
    drawnow;
    DATA = cmb.ReClassifyAll(DATA);
    toc;
end

%   DATA.AllData.Spikes = Spikes;  %% Don't change main probe array
set(DATA.toplevel,'UserData',DATA);
cmb.AddProbeList(DATA);
if isfield(DATA,'currentexpt') & ismember(DATA.probe, probelist)
    id = find(probelist == DATA.probe);
    if id < length(probelist)
        id = probelist(id+1);
    elseif id > 1
        id = probelist(id-1);
    else
        id = DATA.probe;
    end
    xc = cmb.CalcXcorr(DATA,DATA.currentexpt(1),DATA.probe,id);
end

