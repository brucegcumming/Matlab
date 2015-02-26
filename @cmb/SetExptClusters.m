function savecl = SetExptClusters(caller,b, varargin)
%DATA = cmb.combine('getstate');

savecl = 1;
DATA = GetDataFromFig(caller);
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'nosave',5)
        if length(DATA.probes) > 2
            savecl = 0;
        end
    end
    j = j+1;
end

nc = cmb.CountClusters(DATA.cluster);
DATA.nclusters = nc;
eid = DATA.currentexpt(1);
DATA.Expts{eid}.Cluster{1,DATA.probe}.touched = 2;
DATA.Expts{DATA.currentexpt(1)}.gui.ncluster = nc;
DATA.cluster{DATA.currentcluster,DATA.probe}.autocut = 0;
if isfield(DATA.cluster{DATA.currentcluster,DATA.probe},'deleted') & ...
        DATA.cluster{DATA.currentcluster,DATA.probe}.deleted == 1
    DATA.cluster{DATA.currentcluster,DATA.probe}.quality = 0;
else
    DATA.cluster{DATA.currentcluster,DATA.probe}.quality = DATA.state.currentclusterquality;
end

touched = zeros(size(DATA.probelist));
if isfield(DATA,'AllClusters')
    for k = 1:nc
        for j = 1:min([size(DATA.cluster,2) length(DATA.AllClusters)])
            if  isfield(DATA.cluster{k,j}, 'touched') && DATA.cluster{k,j}.touched > 0 && isfield(DATA.AllSpikes{j},'cx')
                DATA.Expts{DATA.currentexpt(1)}.Cluster{k,j} = DATA.cluster{k,j};
                if sum(ismember(DATA.cluster{k,j}.params,[29 30 31])) %templates used
                    DATA.Expts{DATA.currentexpt(1)}.Cluster{k,j}.Templates = DATA.Templates;
                end
                touched(j) = 1;
            end
        end
    end
    fprintf('Probes %s\n',num2str(find(touched >0)));
elseif isfield(DATA,'AllSpikes') % only save clusters that have been set
    for k = 1:nc
        for j = 1:min([size(DATA.cluster,2) length(DATA.AllSpikes)])
            if  isfield(DATA.cluster{k,j}, 'touched') && DATA.cluster{k,j}.touched > 0 && isfield(DATA.AllSpikes{j},'cx')
                DATA.Expts{DATA.currentexpt(1)}.Cluster{k,j} = DATA.cluster{k,j};
                if sum(ismember(DATA.cluster{k,j}.params,[29 30 31])) %templates used
                    DATA.Expts{DATA.currentexpt(1)}.Cluster{k,j}.Templates = DATA.Templates;
                end
                touched(j) = 1;
            end
        end
    end
    fprintf('Probes %s\n',num2str(find(touched >0)));
else
    DATA.Expts{DATA.currentexpt(1)}.Cluster = DATA.cluster;
    for k = 1:nc
        if isfield(DATA.cluster{k,1},'params') & sum(ismember(DATA.cluster{k,1}.params,[29 30 31])) %templates used
            DATA.Expts{DATA.currentexpt(1)}.Cluster{k,1}.Templates = DATA.Templates;
        end
    end
    %
    %  If the set button is hit, it means that user is happy with cluster
    %  currently shown in the XY plot. Need to track this so that ONLY these
    %  clusters are recorded in OnlineClusters. Set touched to 2, then onlhy
    % Add these to AllClusters below
    DATA.Expts{eid}.Cluster{1,DATA.probe}.touched = 2;
end
DATA.Expts{DATA.currentexpt(1)}.gui.clustertype = 1;
DATA.Expts{DATA.currentexpt(1)}.gui.classified = 1;
spkvarnames= DATA.spkvarnames;
p = DATA.probe;
for j = 1:length(DATA.Expts)
    if isfield(DATA.Expts{j},'Cluster')
        enc =  cmb.CountClusters(DATA.Expts{j}.Cluster);
        if DATA.state.online == 1
            for k = 1:size(DATA.Expts{j}.Cluster,2)
                %             for k = p;
                if ~isfield(DATA.Expts{j}.Cluster{1,k},'touched')
                    DATA.Expts{j}.Cluster{1,k}.touched = 0;
                end
                %online temp cuts "just top see" get save here. User needs to delete diligenlyt
                %tried only saving current probe, but that is no good, esp for Grid Data
                if DATA.Expts{j}.Cluster{1,k}.touched > 0
                    AllClusters{j}.Cluster(1:enc,k) = DATA.Expts{j}.Cluster(1:enc,k);
                    touched(p) = 1;
                end
            end
        else
            AllClusters{j}.Cluster = DATA.Expts{j}.Cluster(1:enc,:);
        end
        AllClusters{j}.ids = [DATA.Expts{j}.Trials(1).id DATA.Expts{j}.Trials(end).id];
        clustertypes(j) = DATA.Expts{j}.gui.clustertype;
    end
    excludelist{j}  = find([DATA.Expts{j}.Trials.Trial] < 0);
    if isfield(DATA.Expts{j}.gui,'spks')
        expispk{j} = DATA.Expts{j}.gui.spks;
    else
        expispk{j} = [];
    end
    ispklen(j) = length(expispk{j});
end

if savecl
    a = findobj(DATA.xyfig,'String','Set+Next');
    set(a,'backgroundcolor','w')
    DATA.savedclusters = 1;
else
    DATA.savedclusters = 0;
end
Templates = DATA.Templates;
if isfield(DATA,'TemplateInfo')
    TemplateInfo = DATA.TemplateInfo;
else
    TemplateInfo = [];
end

if DATA.currentcluster > 7 %defining artifact
    id = cmb.PlotArtifacts(DATA);
    if length(id)
        if isfield(DATA,'AllSpikes')
        elseif isfield(DATA.AllData,'Spikes')
            DATA.AllData.Spikes.codes(id,2) = 8;
        end
    end
end
if isfield(DATA,'AllSpikes')  %need to save all probes with a new cluster
    for j = find(touched > 0)
        if isfield(DATA.AllSpikes{j},'times')
            f = DATA.AllSpikes{j}.firstspki;
            clid(f:f+size(DATA.AllSpikes{j}.times,1)-1) = DATA.AllSpikes{j}.codes(:,2);
            save(cmb.ClusterFile(DATA,'probe',j),'AllClusters','clustertypes','excludelist','clid','spkvarnames','Templates','TemplateInfo');
        end
    end
    if strncmp(DATA.filetype,'Grid',4)
        for j = 1:length(DATA.AllSpikes)
            if isempty(DATA.AllSpikes{j})
                codes{j} = [];
                times{j} = [];
            else
                codes{j} = DATA.AllSpikes{j}.codes(:,2);
                times{j} = DATA.AllSpikes{j}.times;
            end
        end
        save(cmb.ClusterFile(DATA,'codes'),'codes','times');
    end
elseif isfield(DATA,'AllClusters')  %need to save all probes with a new cluster
    for j = find(touched > 0)
        clid = DATA.AllClusters(j).codes(:,2);
        save(cmb.ClusterFile(DATA),'AllClusters','clustertypes','excludelist','clid','spkvarnames','Templates','TemplateInfo');
    end
elseif isfield(DATA.AllData.Spikes, 'codes')
    if savecl
        clid = DATA.AllData.Spikes.codes(:,2);
        cfile = cmb.ClusterFile(DATA);
        save(cfile,'AllClusters','clustertypes','excludelist','clid','spkvarnames','ispklen','Templates','TemplateInfo');
        ifile = regexprep(cfile,'\.p([0-9]*)cl\.','.p$1ispk.');
        if ~exist(ifile)
            save(ifile,'expispk','ispklen')
        end
    end
else
    save(cmb.ClusterFile(DATA),'AllClusters','clustertypes','excludelist','clid','spkvarnames','ispklen','Templates','TemplateInfo');
end

if length(DATA.probelist) > 1 && savecl > 1
    save(cmb.ClusterFile(DATA,'allprobes'),'AllClusters','spkvarnames','Templates','TemplateInfo');
end

if DATA.logfid
    fprintf(DATA.logfid,'Cluster set for P%d Expt %d %s by %s\n',DATA.probe,DATA.currentexpt(1),datestr(now),DATA.user);
end
% now re-do list of su-expts to reflect cut clusters
eid = get(DATA.clst,'value');
if DATA.currentcluster > 7 %defining artifact
    DATA.currentcluster = 1;
    figure(DATA.xyfig);
    hold off;
    DATA = cmb.DrawXYPlot(DATA, DATA.Expts{DATA.currentexpt(1)}.gui.spks);
    cid = findobj('Tag','Clusterid');
    set(cid,'value',1);
end
DATA = cmb.ListSubExpts(DATA,eid,'relist');
set(DATA.toplevel,'UserData',DATA);
cid = findobj('Tag','ClusterIsSet');
set(cid,'value',1);
if DATA.state.autoplotcells
    GetFigure(DATA.tag.celllist);
    cmb.PlotCellList(DATA);
end


