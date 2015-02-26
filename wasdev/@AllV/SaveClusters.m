function [DATA, id] = SaveClusters(DATA, outname,varargin)        if nargin == 1        outname = AllV.ClusterFile(DATA.name,DATA.Expt,'subdir',DATA.clustersubdir);    end    quickmode = 0;   if DATA.loadfromspikes    savexy = 2;   else    savexy = 1;   end    j = 1;    while j <= length(varargin)        if strncmpi(varargin{j},'quick',5)            quickmode = 1;            if DATA.quicksave.QuickDetails  == 0                savexy = 0;            end        end        j = j+1;    end        startts = now;    DataClusters = AllV.mygetappdata(DATA,'Clusters');    dname = strrep(outname,'.mat','Details.mat');    logname = [DATA.datadir '/' 'ClusterLogExpt' num2str(DATA.exptno) '.mat'];    if DATA.toplevel    set(DATA.toplevel,'name',sprintf('Saving %s',outname));    drawnow;    end    Vall = AllV.mygetappdata(DATA,'Vall');    if ~isfield(DATA.cluster,'next') || ~iscell(DATA.cluster.next) %should not need this - check        DATA = AllV.AddErr(DATA,'ERROR!!!!!!!!!!! - missing .next\n');        DATA.cluster.next  = {};    end    if quickmode == 0 && ~strcmp(DATA.cluster.automode,'james')          DATA = AllV.ClassifyAll(DATA,0); %Check all are up to date          if DATA.toplevel; set(0,'currentfigure',DATA.toplevel); end    end        if sum(strcmp(DATA.DataType,{'Spike2'  'Default'})) && isempty(DATA.ArrayConfig)        SaveArrayConfig(DATA);    end    DATA = AllV.CheckClusterMarks(DataClusters,DATA);    if DATA.checkclusters    AllV.CheckClusters(DataClusters,'Save');    AllV.CheckClusters(DataClusters,'CheckNexts','Save');    end%Saving should only change one probe. So read in file from disk, and%only update the current probe, then write out.     if exist(outname,'file')        load(outname);        ts = now;        for j = 1:length(DataClusters);            if j > length(Clusters) || isempty(Clusters{j})                Clusters{j} = DataClusters{j};            elseif ~isfield(Clusters{j},'mahal') %empty really                Clusters{j} = DataClusters{j};            elseif isfield(Clusters{j},'mean') %old                Clusters{j} = rmfield(Clusters{j},'mean');            elseif isfield(Clusters{j},'r')                 Clusters{j} = rmfield(Clusters{j},'r');            elseif isfield(Clusters{j},'clst')                 Clusters{j} = rmfield(Clusters{j},'clst');            elseif isfield(Clusters{j},'t')                 Clusters{j} = rmfield(Clusters{j},'t');            end                        if isfield(Clusters{j},'next') && ~iscell(Clusters{j}.next) %old style                last = Clusters{j}.next;                Clusters{j} = rmfield(Clusters{j},'next');                Clusters{j}.next{1} = rmfields(last,'next'); %get rid of next.next            elseif ~isfield(Clusters{j},'next') && isfield(Clusters{j},'mahal')                Clusters{j}.next = {};            end                                    if ~isfield(Clusters{j},'probe')                Clusters{j}.probe = j;            end            if ~isfield(Clusters{j},'auto')                if isfield(DataClusters{j},'auto')                    Clusters{j}.auto = DataClusters{j}.auto;                else                    Clusters{j}.auto = 0;                end                                elseif strfind(outname,'AutoClusterTimes')                Clusters{j}.auto = 1;            end            if ~isfield(Clusters{j},'mahal') && isappdata(DATA.toplevel,'AutoClusters') %still empty                AutoClusters = getappdata(DATA.toplevel,'AutoClusters');                if length(AutoClusters) >= j                    Clusters{j} = AutoClusters{j};                    fprintf('Cluster %d was empty - reverting to Auto\n',j);                    errordlg(sprintf('Cluster %d was empty Reloaded AutoCluster\n',j),'Cluster Error','modal');                                    elseif DATA.checkclusters                    errordlg(sprintf('Cluster %d was empty and no AutoCluster\n',j),'Cluster Error','modal');                                    end                                end        end    else        Clusters = DataClusters;    end        if DATA.profiling        fprintf('ReLoaded Clusters at %.3f\n',mytoc(startts));    end            if length(Clusters) >= DATA.probe(1) && isfield(Clusters{DATA.probe(1)},'auto')    wasauto = Clusters{DATA.probe(1)}.auto;    else        wasauto = 1;    end            if savexy == 2         [DATA, ClusterDetails] = AllV.LoadClusterDetails(DATA,'orauto');    elseif exist(dname,'file') && quickmode == 0        load(dname);    elseif isappdata(DATA.toplevel,'ClusterDetails') && DATA.interactive >= 0        ClusterDetails = AllV.mygetappdata(DATA,'ClusterDetails');        for j = 1:length(ClusterDetails)            if isempty(ClusterDetails{j})                mycprintf('errors','Missing ClsuterDetails %d - will not save\n',j);                savexy = 0;            end        end    elseif savexy > 0 && exist(dname,'file') && exist('ClusterDetails','var') %don't overwrite if missing some cells        savexy = 0;    end    if savexy > 0 && ~exist('ClusterDetails','var')        [DATA, ClusterDetails] = AllV.LoadClusterDetails(DATA,'orauto');    end    if DATA.profiling        fprintf('Got Details at %.3f\n',mytoc(startts));    end        p = AllV.ProbeNumber(DATA);    if DATA.autorefine > 0        DATA.cluster.manual = 3;        res.refinemode = DATA.refinemode;    elseif DATA.cluster.auto ==1 || DATA.cluster.auto == 2        DATA.cluster.manual = 0;    else        DATA.cluster.manual = 1;    end% Check DATA.cluster.next before copying to Clusters{p}    for j = 1:length(DATA.cluster.next)        DATA.cluster.next{j} = rmfields(DATA.cluster.next{j},'r','clst','rescaled');        %Keep empty nexts empty        if isfield(DATA.cluster.next{j},'space')            DATA.cluster.next{j}.manual = DATA.cluster.manual;        end        AllV.CheckClusterValues(DATA, DATA.cluster.next{j});    end      %if this is the first lap of an autocluster we don't want to keep the stored clsuter.first -%what happend last time around    if isempty(DATA.lastcut) && length(Clusters) >= p && isfield(Clusters{p},'first')        Clusters{p} = rmfield(Clusters{p},'first');    end    Clusters{p}.forceevec = 0; %%in case not a field in DATA.cluster    f = fields(DATA.cluster);    for j = 1:length(f)        Clusters{p}.(f{j}) = DATA.cluster.(f{j});    end    Clusters{p}.probe = p;        xid = {};    if isfield(DATA,'clst')        id = find(DATA.clst(DATA.uid) == 2);        nc = length(unique(DATA.clst(DATA.uid)));        for j = 3:nc        xid{j-2} = find(DATA.clst(DATA.uid) == j);        end    else        id = 1:length(DATA.uid);    end    [a,b] = Counts(DATA.clst(DATA.uid));    if length(a) < 2        a(2) = 0;    end    DATA.savespkid = DATA.uid;    fprintf('Saving (%s)/%d Spikes (P%d) to %s (%.2f)\n',sprintf('%d ', a(2:end)),length(DATA.uid),p,outname,mytoc(startts));% ClusterDetails records all event times, and the classification (clst)%Clusters just has the times of the classified  events = smallest file ffor%combine    Clusters{p}.times = DATA.t(DATA.uid(id));    for c = 1:length(Clusters{p}.next)        if isfield(Clusters{p}.next{c},'space')            id = find(DATA.clst(DATA.uid) == 2+c);            Clusters{p}.next{c}.times = DATA.t(DATA.uid(id));        end    end    Clusters{p}.errs = DATA.errs;    Clusters{p}.user = DATA.user;    Clusters{p}.hostname = DATA.hostname;    if quickmode == 0    Clusters{p}.dpsum = sum(abs(DATA.cluster.MeanSpike.vdprime(DATA.chspk,:)),2);    end    if ~isempty(DATA.restricttimerange)        Clusters{p}.restricttimerange = DATA.restricttimerange;    elseif isfield(Clusters{p},'restricttimerange')        Clusters{p} = rmfield(Clusters{p},'restricttimerange');    end    if ~isempty(DATA.excludetrialids)        Clusters{p}.excludetrialids = DATA.excludetrialids;    elseif isfield(Clusters{p},'excludetrialids')        Clusters{p} = rmfield(Clusters{p},'excludetrialids');    end    Clusters{p}.spkfile = AllV.SpkFileName(DATA);    Clusters{p}.exptno = DATA.Expt.exptno;    ClusterDetails{p}.Evec = DATA.Evec;    if isfield(DATA,'rV')        ClusterDetails{p}.triggerV = DATA.rV;    end%when calling reclassify, ClusterDetails is loaded, so can get evec from%there. But its cheap (ish- why not have in in Clusters? %    Clusters{p}.Evec = DATA.Evec;    if isfield(DATA.Expt.Header,'ReadMethod')        DATA.cluster.exptreadmethod = DATA.Expt.Header.ReadMethod;    else        DATA.cluster.exptreadmethod = 0;    end%remove fields that we don't want saved. N.B. might be inherited from an old Clusters file. %xtimes is no longer used%xy, goes into ClusterDetails%Evec more than doubles size, so keep this in ClusterDetails. Only need it%for Reclassify.    Clusters{p} = rmfields(Clusters{p},'r','xtimes', 'rescaled', 'h', 'dragfcn' ,'handles', 'plotargs', 'done' ,'selected' ,'down', 'color', 'axis' );    if Clusters{p}.space(1) == 1        Clusters{p}.Evec = DATA.Evec;        Clusters{p}.Evec.Evec = DATA.Evec.Evec(:,1:5);    else        Clusters{p} = rmfields(Clusters{p},'Evec');    end    Clusters{p}.trigdt = DATA.trigdt;    Clusters{p}.tsmooth = DATA.triggersmooth;    Clusters{p}.triggerchan = DATA.triggerchan;    Clusters{p}.triggertype = DATA.triggertype;    Clusters{p}.clst  = DATA.clst(DATA.uid);    Clusters{p}.version = DATA.version;    Clusters{p}.duration = DATA.duration;    Clusters{p}.probe = p;    Clusters{p}.usealltrials = DATA.usealltrials;        if isfield(DATA,'DataType')        Clusters{p}.DataType = DATA.DataType;    end    if ~isfield(DATA.cluster,'auto')        DATA.cluster.auto = 0;    end    if DATA.userefcluster         Clusters{p}.auto = DATA.cluster.auto;        Clusters{p}.recluster =  3;    elseif DATA.recluster ~= 2        Clusters{p}.auto = 0;        Clusters{p}.recluster = DATA.recluster + 100 * DATA.forcecluster;    else        Clusters{p}.auto = DATA.cluster.auto;        Clusters{p}.recluster =  100 * DATA.forcecluster;    end    Clusters{p}.pcmean = mean(DATA.pcs(:,1:4)); %to check for sign reversal    Clusters{p}.memsz = [AllV.memsize(DATA) DATA.fullvsize];    if isfield(DATA.cluster,'DprimeUsed')        Clusters{p}.DprimeUsed = DATA.cluster.DprimeUsed;    end    if isfield(DATA.cluster,'mumeanUsed')        Clusters{p}.mumeanUsed = DATA.cluster.mumeanUsed;    end    if isfield(DATA.cluster,'TemplateUsed')        Clusters{p}.TemplateUsed = DATA.cluster.TemplateUsed;    elseif isfield(DATA,'TemplateUsed')  && ~isfield(Clusters{p},'TemplateUsed')        Clusters{p}.TemplateUsed = DATA.TemplateUsed;        if isfield(DATA,'DprimeUsed') % can happne if templates not calculated for this probe            Clusters{p}.DprimeUsed = DATA.DprimeUsed;        end    end    Clusters{p}.chspk = DATA.chspk;    Clusters{p}.vsmps = DATA.vsmps;    if ~isfield(Clusters{p},'errs')        Clusters{p}.errs = {};    end    if isfield(DATA.cluster,'errs')        for j = 1:length(DATA.cluster.errs)            id = strcmp(DATA.cluster.errs{j},Clusters{p}.errs);            if isempty(id)                Clusters{p}.errs = {Clusters{p}.errs{:} DATA.cluster.errs{j}};            end        end    end    if isfield(DATA.cluster,'bestspace')    if DATA.cluster.bestspace(2) == 3        Clusters{p}.vspace = DATA.vspace;    end    end    if isfield(DATA,'jamescluster') && strcmp(DATA.autocutmode,'james')        Clusters{p}.jamescluster = DATA.jamescluster;    end    if DATA.checkclusters        AllV.CheckClusters(Clusters,'CheckFitSpace');    end    Clusters{p}.isicheck = DATA.isicheck;    if length(DATA.artifacttimes)        Clusters{p}.artifacttimes = DATA.artifacttimes;    end    Clusters{p}.clusterprog = sprintf('AllVPcs %.2f',DATA.version);    Clusters{p}.progversion = DATA.version;    Clusters{p}.missingtrials = DATA.missedtrials;    if DATA.autorefine > 0         Clusters{p}.manual = 3;    elseif DATA.cluster.auto ==1        Clusters{p}.manual = 0;    else        Clusters{p}.manual = 1;    end    AllV.CheckClusterValues(DATA, DATA.cluster);        Clusters{p} = rmfields(Clusters{p},{'clst' 'r'});    Clusters{p}.addmean = DATA.addmean;    Clusters{p}.savetime(1) = now;    if AllV.NeedTemplateForCluster(Clusters{p},1) ==2 && isfield(DATA.Template,'othermeans')        Clusters{p}.MeanSpike.othermeans = DATA.Template.othermeans(2:end);    end    if DATA.trigdt == 4 && isfield(DATA,'TriggerTemplate')        Clusters{p}.TriggerTemplate = DATA.TriggerTemplate;    end    if isfield(DATA,'xy')        fprintf('%s Setting xy, t, and clst in ClusterDetails (%.2f)\n',AllV.IDStr(DATA),mytoc(startts));            ClusterDetails{p}.xy = DATA.xy{1}(DATA.uid,:);            ClusterDetails{p}.t = DATA.t(DATA.uid);            ClusterDetails{p}.clst = DATA.clst(DATA.uid);        if DATA.interactive >= 0 && isappdata(DATA.toplevel,'ClusterDetails')            DATA = AllV.mysetappdata(DATA,'ClusterDetails',ClusterDetails);            if DATA.auto.showxysaved                AllV.PlotAllProbes(DATA, 'allxy', 'probes',p,'linewidth',2);            end        end    end    ClusterDetails{p}.angle = Clusters{p}.angle;    ClusterDetails{p}.shape = Clusters{p}.shape;    if isfield(Clusters{p},'crit')        ClusterDetails{p}.crit = Clusters{p}.crit;    else        ClusterDetails{p}.crit = NaN;    end            meanvfile = [DATA.datadir '/Expt' num2str(DATA.exptno) 'meanV.mat'];    if ~exist(meanvfile) && isfield(Vall,'meanV');        FullV = rmfield(Vall,'V');        fprintf('Saving %s\n',meanvfile);        save(meanvfile,'-v7.3','FullV');    end    if savexy && isfield(DATA,'xy')        if DATA.savetrigger %triggeV always saved now. Shound not need this            ClusterDetails{p}.rV = DATA.rV;        end        if savexy ==2 || savexy ==1 %used to be only for ==1, but suresly this is wrong            for j = 1:length(Clusters{p}.next)                if length(DATA.xy) > j && ~isempty(Clusters{p}.next{j}) ...                        && ~isempty(DATA.xy{j+1})                    ClusterDetails{p}.next{j}.xy = DATA.xy{j+1}(DATA.uid,:);                end            end        end        id = DATA.uid;        if diff(size(ClusterDetails{p}.clst)) > 0            fprintf('Clst size 2 is %d\n',size(ClusterDetails{p}.clst,2));            ClusterDetails{p}.clst = ClusterDetails{p}.clst';        end        ClusterDetails{p}.ctime = Clusters{p}.ctime;    else        id = [];    end    %if making an automatic cut, and a manual one was made already, dont     %overwrite the spikes file.  To force overwrite, set savespikes to 2.    if DATA.savespikes ==1 && wasauto == 0 && DATA.cluster.auto == 1        fprintf('Will not save spikes(auto) - manual cut was made\n');        id = [];    end    if DATA.profiling        fprintf('Set Details at %.3f\n',mytoc(startts));    end        if quickmode || DATA.savespikes == 3        DataClusters{p} = Clusters{p};    end    Clusters{p}.savetime(2) = now;        if DATA.nolog == 0 && quickmode >= 0 %save or not in quick?        ClusterLog = {};        if exist(logname,'file')            load(logname);            ncl = length(ClusterLog)+1;        else            ncl = 1;        end        ClusterLog{ncl}.cluster = 1;        C = Clusters{p};        ClusterLog{ncl}.shape = C.shape;        ClusterLog{ncl}.space = C.space;        ClusterLog{ncl}.mahal = C.mahal;        ClusterLog{ncl}.quick = quickmode;        if isfield(C,'crit')            ClusterLog{ncl}.crit = C.crit;        else            ClusterLog{ncl}.crit = NaN;        end        ClusterLog{ncl}.angle = C.angle;        if C.shape == 0            ClusterLog{ncl}.xyr = C.xyr;        end        ClusterLog{ncl}.savetime = now;        ClusterLog{ncl}.user = DATA.user;        ClusterLog{ncl}.probe = p;        ClusterLog{ncl}.exptno = DATA.exptno;        ClusterLog{ncl}.hostname = gethostname;        ClusterLog{ncl}.recluster = [DATA.recluster DATA.cluster.auto];        ClusterLog{ncl}.ncut = [C.ncut DATA.nevents];        if isfield(C,'reclassify')            ClusterLog{ncl}.reclassify = C.reclassify;        end        for j = 1:length(C.next)            next.cluster = j+1;            if ~isempty(C.next{j})                ClusterLog{ncl}.next{j} = CopyFields(next,C.next{j},{'xyr', 'angle', 'shape',' space', 'mahal' ,'crit'});            end        end        save(logname,'ClusterLog');        if isfield(DATA.cluster,'gmfit')            nd = size(DATA.cluster.gmfit.mu,2);        else            nd = 0;        end        if DATA.logfid > 2 && DATA.logfid < 20            fprintf(DATA.logfid,'E%dP%d Space%s Fit %d dims Saved (%s)/%d spikes at %s\n',...                DATA.exptno,p,sprintf(' %d',DATA.cluster.space),nd,datestr(now),...                sprintf('%d ',Counts(DATA.clst)),DATA.nevents);        end    end        nerr = 0;    badc = 0;    for j = 1:length(Clusters)        if ~isfield(Clusters{j},'mahal')            nerr = nerr+1;            badc = j;        end    end    if nerr > 0 && DATA.checkclusters        errordlg(sprintf('Cluster %d was empty and no AutoCluster\n',badc),'Cluster Error','modal');    end    if isempty(Vall)        if isfield(DATA,'FullVData')            FullVData = DATA.FullVData;        else            FullVData = [];        end    else        FullVData = MakeFullVInfo(Vall);        DATA.FullVData = FullVData;    end    Clusters{p} = CopyFields(Clusters{p},FullVData,{'coilnoiseratio','highpass'});    FullVData = rmfields(FullVData,'meanV'); %just in case - lots of space    if savexy > 0 && DATA.savespikes ~= 3        if DATA.profiling            fprintf('Saving Details at %.3f\n',mytoc(startts));        end        if strcmp(DATA.autocutmode,'james')            ClusterDetails{p} = CondenseDetails(ClusterDetails{p});        end        ClusterDetails{p}.savetime = now;        try        save(dname,'-v7.3','ClusterDetails','FullVData');        catch            cprintf('Error Saving %s.  Trying delete first\n',dname);            delete(dname);            save(dname,'-v7.3','ClusterDetails','FullVData');        end    else        DATA = AllV.mysetappdata(DATA,'ClusterDetails',ClusterDetails);    end    Clusters{p}.savetime(3) = now;    if DATA.profiling        fprintf('Saved Details at %.3f\n',mytoc(startts));    end        if DATA.auto.uselastcluster    %Keep track of last saved cluster for each probe        for j = 1:length(Clusters)            if j > length(DATA.LastClusters) || isempty(DATA.LastClusters{j})                DATA.LastClusters{j} = Clusters{j};            end        end        DATA.LastClusters{p} = Clusters{p};    end    if DATA.savespikes ==3        DATA = AllV.mysetappdata(DATA,'Clusters',Clusters);        return;    end    if DATA.auto.backupcluster && exist(outname)        [a,b,c] = fileparts(outname);        backdir = [a '/backup'];        if ~exist(backdir,'dir')            mkdir(backdir);            fileattrib(backdir,'+w','g');        end        backfile = [backdir '/' b datestr(now,'mmddyy.HHMMSS') 'P' num2str(p) c];        try            ts = now;            movefile(outname,backfile);            fprintf('Backup took %.2f (%.2f)\n',mytoc(ts),mytoc(startts));        catch ME            DATA = AllV.AddErr(DATA,ME,'-show','Cant move %s to %s - check permissions\n',outname,backfile);        end        mvd = dir(outname);        if ~isempty(mvd)            DATA = AllV.AddErr(DATA,'-show','%s Still exists!!!! Backup should have moved\n',outname);        end    end    if DATA.auto.backupcluster == 0        cprintf('red','Clusters not being backup up before save - NOT RECOMENDED\n');    end        try        save(outname,'-v7.3','Clusters','FullVData');    catch        DATA = AllV.AddErr(DATA,'-show','Error writing %s. Trying delete first\n',outname);        save(outname,'-v7.3','Clusters','FullVData');    end    if DATA.profiling        fprintf('Saved at %.3f\n',mytoc(startts));    end    DATA.savespkid = id;    if DATA.interactive < 0        checksaves = 1;    else        checksaves = 0; %testing erors saving ClusterDetails. Savinging wiht 7.3 about shoule make this unnecessary    end    if savexy > 0 && checksaves;        try        load(dname);        catch            fprintf('Saving %s was corrupt. Trying again\n',dname);            delete(dname);            save(dname,'-v7.3','ClusterDetails','FullVData');           load(dname);        end    end        if DATA.toplevel && DATA.interactive >= 0        DataClusters{p} = Clusters{p};        DATA = AllV.mysetappdata(DATA,'Clusters',DataClusters);        DATA = AllV.mysetappdata(DATA,'ClusterDetails',ClusterDetails);        set(DATA.toplevel,'name',get(DATA.toplevel,'Tag'));    end    