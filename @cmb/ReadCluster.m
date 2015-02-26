function [DATA, D] = ReadCluster(DATA, eid, pid, varargin)

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
ts = now;
if eid > DATA.appendexptid
    [Clusters, F, D] = LoadCluster(DATA.appenddir{1},eid-DATA.appendexptid(1)+1,'getxy','noauto');
else
    [Clusters, F, D] = LoadCluster(DATA.datadir,exptno,'getxy');
end
%NB loadtime is the time at which load was made, used for checking whether
%file has changed. Loaddur is time taken to do load
D.loaddur = mytoc(ts);
if isfield(D,'err') && ~isempty(D.err)
    acknowledge(D.err,DATA.toplevel);
end
lasterr = '';
cid = 0;
if isfield(DATA.CellDetails,'trialids')
    if checkcl
        try
            DATA.CellDetails = CheckExptClusters(Clusters,DATA.Expts,DATA);
        catch ME
            CheckExceptions(ME);
        end
    end
    cid = find(DATA.CellDetails.exptids == eid); %cid is index to expt list in Celllist.
elseif checkcl
    cprintf('blue','Celldetails has no trialids\n');
end

diperrs = {};
for j = length(Clusters):-1:1
    
    if isfield(Clusters{j},'xy')
        if isfield(Clusters{j},'trigset')
            if isfield(Clusters{j}.trigset{1},'t')                
                allt = [Clusters{j}.t Clusters{j}.trigset{1}.t];
                allcl = [Clusters{j}.clst; Clusters{j}.trigset{1}.clst];
            else
                allt = [Clusters{j}.t];
                allcl = [Clusters{j}.clst];
            end
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
        if isfield(Clusters{j},'xyr')
            DATA.AllClusters{eid}(j).xyr = Clusters{j}.xyr;
            DATA.AllClusters{eid}(j).angle = Clusters{j}.angle;
            if isfield(Clusters{j},'aspectratio')
                DATA.AllClusters{eid}(j).aspectratio = Clusters{j}.aspectratio;
            end
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
                    elseif regexp(err,'Not Applying RefCluster')
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
                if isfield(Clusters{j}.next{k},'xy')
                    DATA.AllClusters{eid}(j).next{k}.xy = Clusters{j}.next{k}.xy;
                else
                    cprintf('red','E%dP%d cluster %d No xy data - using cluster 1\n',eid,j,k+1);
                    DATA.AllClusters{eid}(j).next{k}.xy = Clusters{j}.xy;
                end
                    
                if isfield(Clusters{j}.next{k},'xyr')
                    DATA.AllClusters{eid}(j).next{k}.xyr = Clusters{j}.next{k}.xyr;
                    DATA.AllClusters{eid}(j).next{k}.angle = Clusters{j}.next{k}.angle;
                    if isfield(Clusters{j}.next{k},'xyr')
                        DATA.AllClusters{eid}(j).next{k}.aspectratio = Clusters{j}.next{k}.aspectratio;
                    end
                else
                    fprintf('No xy data \n');
                end
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
        DATA.AllClusters{eid}(j).codes = [];
    end
end
if ~isfield(DATA,'Clusterfile') || length(DATA.Clusterfile) < eid
    DATA.Clusterfile{eid} = D;
else
    DATA.Clusterfile{eid} = CopyFields(DATA.Clusterfile{eid},D);
end
DATA.Expts{eid}.gui.counted = 0;

