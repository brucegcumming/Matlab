function [DATA, D] = LoadCluster(DATA, eid, pid, varargin)

checkcl = 0; 
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'check',5)
        checkcl =1;
    end
    j = j+1;
end

if isfield(DATA,'suffixlist')
    exptno = DATA.suffixlist(eid); %suffix #
else %Grid files
    exptno = eid;
end
if eid > DATA.appendexptid
    [Clusters, F, D] = LoadCluster(DATA.appenddir{1},eid-DATA.appendexptid(1)+1,'getxy','noauto');
else
    [Clusters, F, D] = LoadCluster(DATA.datadir,exptno,'getxy');
end
lasterr = '';
cid = 0;
if isfield(DATA.CellDetails,'trialids')  && checkcl
    try
        DATA.CellDetails = CheckExptClusters(Clusters,DATA.Expts,DATA);
    catch ME
        CheckExceptions(ME);
    end
    cid = find(DATA.CellDetails.exptids == eid);
    if ~isempty(cid)
        cid = DATA.CellDetails.exptids(cid);
    end
else
    cprintf('blue','Celldetails has no trialids\n');
end

for j = length(Clusters):-1:1
    
    if isfield(Clusters{j},'xy')
        if isfield(Clusters{j},'trigset')
            allt = [Clusters{j}.t Clusters{j}.trigset{1}.t];
            allcl = [Clusters{j}.clst; Clusters{j}.trigset{1}.clst];
            clist = ClusterDef(Clusters{j}.trigset{1});
            clist = clist(clist>1);
            clist = clist(1)-1;
            allxy = cat(1,Clusters{j}.xy,Clusters{j}.next{clist}.xy);
            [t,tid] = sort(allt);
            DATA.AllClusters{eid}(j).times = t' .*10000;
            DATA.AllClusters{eid}(j).cx = allxy(tid,1);
            DATA.AllClusters{eid}(j).cy = allxy(tid,2);
            DATA.AllClusters{eid}(j).codes= allcl(tid)-1;
        else
            DATA.AllClusters{eid}(j).cx = Clusters{j}.xy(:,1);
            DATA.AllClusters{eid}(j).cy = Clusters{j}.xy(:,2);
            DATA.AllClusters{eid}(j).times = Clusters{j}.t' .*10000;
            DATA.AllClusters{eid}(j).codes= Clusters{j}.clst-1;
        end
        DATA.AllClusters{eid}(j).suffix = exptno;
        [DATA.AllClusters{eid}(j).dips, diperrs{j,1}] = cmb.GetDips(Clusters{j});
        DATA.AllClusters{eid}(j).dropi = Clusters{j}.dropi(3);
        if isfield(Clusters{j},'excludetrialids')
            DATA.AllClusters{eid}(j).excludetrialids{1} = Clusters{j}.excludetrialids;
        else
            DATA.AllClusters{eid}(j).excludetrialids{1} = [];
        end
        if isfield(Clusters{j},'missingtrials')
            DATA.AllClusters{eid}(j).missingtrials = Clusters{j}.missingtrials;
        else
            DATA.AllClusters{eid}(j).missingtrials = [];
        end
        if isfield(Clusters{j},'errs')
            showerror = [];
            for e = 1:length(Clusters{j}.errs)
                err = Clusters{j}.errs{e};
                idstr = sprintf('E%dP%d',Clusters{j}.exptno,Clusters{j}.probe(1));
                if ~strcmp(err,lasterr)
                    if regexp(err,'Removed [0-9]* suspic')
                    elseif err(end) == 10
                        fprintf('%s%s',idstr,err);
                    else
                        fprintf('%s%s\n',idstr,err);
                    end
                    if strncmp(err,'Missing Trials',14)
                        showerror = [showerror err];
                    end
                    lasterr = err;
                end
            end
        end
    else
        DATA.AllClusters{eid}(j).cx = [];
        DATA.AllClusters{eid}(j).cy = [];
        DATA.AllClusters{eid}(j).times = [];
        DATA.AllClusters{eid}(j).dips = [];
        DATA.AllClusters{eid}(j).dropi = [];
        DATA.AllClusters{eid}(j).excludetrialids{1} = [];
    end
    if isfield(Clusters{j},'user')
        DATA.AllClusters{eid}(j).user = Clusters{j}.user;
    end
    if isfield(Clusters{j},'savetime')
        DATA.AllClusters{eid}(j).savetime = Clusters{j}.savetime(end);
    end
    
    if isfield(Clusters{j},'next')
        for k = 1:length(Clusters{j}.next)
            if isfield(Clusters{j}.next{k},'dropi') %don't read empties
                [DATA.AllClusters{eid}(j).next{k}.dips, diperrs{j,k+1}] = cmb.GetDips(Clusters{j}.next{k});
                DATA.AllClusters{eid}(j).next{k}.dropi = Clusters{j}.next{k}.dropi(3);
            end
        end
    end
    if isfield(Clusters{j},'mahal')
        DATA.Expts{eid}.Cluster{j}.mahal = Clusters{j}.mahal;
    end
    if cid
        xcl = DATA.CellDetails.excludetrials(cid,j,:);
        oldx = DATA.AllClusters{eid}(j).excludetrialids;
        for k = 1:length(xcl)
            if ~isempty(xcl{k})
                if k > length(oldx)
                    oldx{k} = [];
                end
                DATA.AllClusters{eid}(j).excludetrialids{k} = [oldx{k} xcl{k}];
            end
        end
    end
end
dipmissing = find(CellToMat(diperrs,'probe') > 0);
if ~isempty(dipmissing)
    fprintf('%s Probes %s',diperrs{dipmissing(1)}.s,sprintf('%d,',dipmissing));
end

if isempty(Clusters)
    for j = length(DATA.probelist):-1:1
        DATA.AllClusters{eid}(j).cx = [];
        DATA.AllClusters{eid}(j).cy = [];
        DATA.AllClusters{eid}(j).times = [];
        DATA.AllClusters{eid}(j).dips = [];
    end
end
if ~isfield(DATA,'Clusterfile') || length(DATA.Clusterfile) < eid
    DATA.Clusterfile{eid} = D;
else
    DATA.Clusterfile{eid} = CopyFields(DATA.Clusterfile{eid},D);
end
DATA.Expts{eid}.gui.counted = 0;

