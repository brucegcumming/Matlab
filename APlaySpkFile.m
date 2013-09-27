function [Expt, Expts, AllData, Raw] = APlaySpkFile(name, varargin)

% [Expt, Expts, AllData] = APlaySpkFile(name, ...)
%
% Builds/retnes lists of trials/expts from a matlab file
% APlaySpkFile(name, 'relist')  Re-builds the index filt
% APlaySpkFile(name, 'setprobe', n)  sets which probes data is loaded
% initially
% APlaySpkFile(name, 'setprobe', -1) rebuilds probe list


SpkDefs;
playspikes = 0;
defaults.fz = 96;
starttrial = 0;
onlinedata = 0;
idxfile = [];
Expts = [];
Expt = [];
AllData = [];
s2version = 0;
testframe = 0;
findtrial = 0;
spkch = 'Ch5';
stimch = 'Ch8';
framechname = 'Ch7';
mainsname = 'Ch24';
dstimch = 'Ch18';
ustimmarkname = '';
ustimmarkch = [];
mainsch = [];
stimlvl = [];
logfid = 0;
fixup = 1;
UstimV = 0;
setprobe = 1;
state.nospikes = 0;
nerr = 0;
preperiod=2000;
postperiod = 2000;
savedvdt = 0;
dvfile = [];
Raw = [];
timeoffset = 0;
ignoreSpikeO = 0;  %% for files with full voltage records elsewhere, ignore these
quickload = 0;
%method for matching up StimChan starts with text events. With many channel
%recordings, can get big delays between these, so the original method,
%finding the nearest stimchan event, is unreliable. New method (1) finds
%all events firts, and if the numbers are the same, matches in order.
%dfeault set to method 1 Jan 6 2012 by bgc
state.method = 1; 
state.showerrs = 1;
saveexpts = 0;
state.needframes = 1;
state.nospikes = 0;
state.alltrials = 0;
mkidx =0;  %%need this up here so taht relist works

j = 1;
while j <= length(varargin)
    vg = varargin{j};
    if strncmpi(vg,'alltrials',4)
        onlinedata = 1;
    elseif ischar(vg) & strncmpi(vg,'online',4)
        onlinedata = 1;
    elseif ischar(vg) & strncmpi(vg,'setprobe',4)
        j = j+1;
        setprobe = varargin{j};
    elseif strncmpi(vg,'findprobes',6)
       mkidx = 1;
       setprobe = -1; %force relisting of probes
    elseif strncmpi(vg,'noerrs',5)
        state.showerrs = 0;
    elseif strncmpi(vg,'nospikes',5)
        state.nospikes = 1;
    elseif strncmpi(vg,'rfs',3)
       state.nospikes = 2;
    elseif strncmpi(vg,'bysuffix',7)
       ignoreSpikeO = 2;
    elseif strncmpi(vg,'method',6)
        j = j+1;
        state.method = varargin{j};
    elseif strncmpi(vg,'noframes',6)
        state.needframes = 0;
    elseif strncmpi(vg,'quicksuffix',9)
       ignoreSpikeO = 2;
       quickload = 1;
    elseif strncmpi(vg,'relist',4)
        mkidx = 1;
    elseif strncmpi(vg,'saveexpts',6)
        saveexpts = 1;
    elseif strncmpi(vg,'sortexpts',6)
        idxfile = strrep(name,'.mat','idx.mat');
        load(idxfile);
        [Expts, Expt] = SortExpts(Expts, Expt.Trials, Expt.Header,1, Expt, state);
        return;
    elseif strncmpi(vg,'usealltrials',8)
        state.alltrials = 1;
    elseif strncmpi(vg,'timeoffset',8)
        j = j+1;
        timeoffset = varargin{j};
    end
    
    j = j+1;
end
    
clusterdate = now;
thecluster = 1;

if ischar(name)
    if ~exist(name,'file')
        fprintf('No file %s\n',name);
        return;
    end
    
    
    if strfind(name,'idx.mat')
        argon = {};
        load(name);
        [Expts, Expt] = SortExpts(Expts, Expt.Trials, Expt.Header, thecluster, Expt, state, argon{:});
        return;
    end
    np = 0;
    nlfp = 0;
    nspkt = 0;
    probes = [];
    logname = strrep(name,'.mat', '.log');
%    fprintf('Log %s\n',logname);
   logfid = fopen(logname,'a');
   Oprobe = 0;
   if onlinedata
       %        oname = strrep(name,'online','online2');
       oname = strrep(name,'/Expt','A/Expt');
       oname = regexprep(oname,'(\.[0-9]*.mat)','A$1');
       idxfile = strrep(name,'.mat','idx.mat');
       if exist(idxfile,'file') && ~mkidx
           load(idxfile);
           probes = Expt.Probes;
       else
       if exist(oname,'file')
           af = load(oname);
           f = fields(af);
       else
           f = {};
       end
       for j = 1:length(f)
           if ~isempty(regexp(f{j},'Ch[0-9]*'))
               ch = af.(f{j});
               if strncmpi(ch.title,'Spike',5)
                   np = np+1;
                   if strncmpi(ch.title,'SpikeO',6)
                       probe = sscanf(ch.title,'SpikeO%d');
                       Oprobe = Oprobe+1;
                   else
                       probe = sscanf(ch.title,'Spike %d');
                   end
                   if isempty(probe)
                       probes(np).probe = sscanf(vars{j},'Ch%d');
                   else
                       probes(np).probe = probe;
                   end
                   probes(np).var = f{j};
                   probes(np).traces = ch.traces;
                   probes(np).source = 2;
                   if probe == setprobe(1)
                       Chspk = ch;
                       spkch = 'Chspk';
                   end
               elseif strncmpi(ch.title,'4Trode',5)
                   np = np+1;
                   probe = sscanf(ch.title,'4Trode%d');
                   if isempty(probe)
                       probes(np).probe = sscanf(vars{j},'Ch%d');
                   else
                       probes(np).probe = probe;
                   end
                   probes(np).var = f{j};
                   probes(np).traces = ch.traces;
                   probes(np).source = 2;
                   if probe == setprobe(1)
                       Chspk = ch;
                       spkch = 'Chspk';
                   end
               end
           end
       end
       end
   end
%    fprintf('Reading %s\n',name);
   mkmatver = 0;

 
   if ignoreSpikeO == NaN % don't need this any more, for online at least
       oname = regexprep(name,'.([0-9]*.mat)','A.$1');
       if exist(oname,'file')
           load(oname);
       end
       avars = who('Ch[0-9]*');
       for j = 1:length(avars)
           eval([avars{j} 'A = ' avars{j} ';']);
           clear(avars{j});
       end

   end

   load(name);
   
    
    if exist('SMRFiles','var') %% Concatentae existing files
        toff = 0;
        for j = 1:length(SMRFiles.Names)
            [a,b,c] = APlaySpkFile(SMRFiles.Names{j},varargin{:},'timeoffset',toff);
            if j == 1
                Expt = a;
                Expts = b;
                AllData = c;
            else
                f = fields(Expt.Trials);
                for k = 1:length(f)
                    if isfield(a.Trials,f{k})
                        if diff(size(Expt.Trials.(f{k}))) > 1
                            dim = 2;
                        else
                            dim = 1;
                        end
                    Expt.Trials.(f{k}) = cat(dim,Expt.Trials.(f{k}), a.Trials.(f{k}));
                    end
                end
                Expt.Probes = [Expt.Probes a.Probes];
                Expt.Spkid = cat(1,Expt.Spkid,a.Spkid);
                Expts = [Expts b];
                f = fields(AllData.Spikes);
                for k = 1:length(f)
                AllData.Spikes.(f{k}) = cat(1,AllData.Spikes.(f{k}),c.Spikes.(f{k}));
                end
            end
            toff = max(Expt.Trials.End)+1;

        end
        return;
    end
    if state.nospikes == 2  %make ufl file
        MkUfl(name,Ch30,'overwrite');
        return;
    end
    if exist('Ch31','var') & isfield(Ch31,'comment')
        mkmatver = sscanf(Ch31.comment,'MkMat V%n');
    end
    if exist('Ch30','var') & strncmp(Ch30.comment,'GridData',8) 
        state.nospikes = 1; %Don't try to load up all spike files if its Utah Array
    end

    vars = who('Ch*');
    for j = 1:length(vars)
        if ~isempty(regexp(vars{j},'Ch[0-9][0-9]*'))
            eval(['ch = ' vars{j} ';']);
            chn = sscanf(vars{j},'Ch%d');
            if chn > 400  %a memory/extra offline channel
            elseif strncmpi(ch.title,'SpikeO',6) && ignoreSpikeO;
            elseif strncmpi(ch.title,'Spike',5) && state.nospikes == 0
                np = np+1;  
                   if strncmpi(ch.title,'SpikeO',6)
                       probe = sscanf(ch.title,'SpikeO%d');
                       Oprobe = Oprobe+1;
                   else
                       probe = sscanf(ch.title,'Spike %d');
                   end
                if isempty(probe)
                    probes(np).probe = sscanf(vars{j},'Ch%d');
                    probes(np).probech = sscanf(vars{j},'Ch%d');
                else
               probes(np).probe = probe;
                end
                probes(np).var = vars{j};
                probes(np).traces = ch.traces;
                probes(np).source = 1;


            elseif strncmpi(ch.title,'4Trode',5)
                np = np+1;  
                probe = sscanf(ch.title,'4Trode%d');
                if isempty(probe)
                    probes(np).probe = sscanf(vars{j},'Ch%d');
                else
               probes(np).probe = probe;
                end
                probes(np).var = vars{j};
                if probe == setprobe(1)
                    Chspk = ch;
                    spkch = 'Chspk';
                end
                probes(np).traces = ch.traces;
                probes(np).source =1;

            elseif strncmpi(ch.title,'uStimMk',7)
                ustimmarkame = vars{j};
                ustimmarkch = ch;
            elseif strncmpi(ch.title,'uStim',5)
                np = np+1;  
                probe = 100;
                UstimV = ch;
                if isempty(probe)
                    probes(np).probe = sscanf(vars{j},'Ch%d');
                    probes(np).probech = sscanf(vars{j},'Ch%d');
                else
               probes(np).probe = probe;
                end
                probes(np).var = vars{j};
                
            elseif strncmpi(ch.title,'StimOn',6)
                stimch = vars{j};
                stimlvl = ch;
            elseif strncmpi(ch.title,'StimChan',8) % stim change detector
                dstimch = vars{j};
                stimchange = ch;
            elseif strncmpi(ch.title,'VTR',3)
                framechname = vars{j};
                framech = ch;
   %             fprintf('Frames in %s\n',vars{j});
            elseif strncmpi(ch.title,'Mains',5)
                mainsname = vars{j};
                mainsch = ch;
            elseif strncmpi(ch.title,'DigMark',7)
                ustimmarkame = vars{j};
                ustimmarkch = ch;
            end
        end
    end
    if np > 0
       if np == 1
           probes(1).probe = 1;
       end
        if exist('Chspk','var') %found selected probe
            spkch = 'Chspk';
        elseif setprobe(1) > 0 
            id = find([probes.probe] == setprobe(1));
            if ~isempty(id)
             spkch = probes(id(1)).var;
            else
             spkch = probes(1).var;
            end
        else
        spkch = probes(1).var;
        end
    else
        idxfile = strrep(name,'.mat','probes.mat');
        if exist(idxfile,'file') & setprobe(1) >= 0
            load(idxfile);
            id = find([probes.probe] == setprobe(1));
            [dp, pref] = fileparts(name);
            spkdir = [dp '/Spikes'];
 % if IgnoreSpikeO is 2 this means don't load the spike waveforms, just the cluster
 % cluster xy values from the ClusterDetails File
            if ~exist(spkdir,'dir') || ignoreSpikeO == 2
                spkdir = dp;
            end

            
            if length(id) > 1 
                ch = sscanf(probes(id(1)).var,'Ch%d');
            else
                ch = sscanf(probes(1).var,'Ch%d');
            end

            if isfield(probes,'filename')
                pname = [spkdir '/' probes(id(1)).filename];
            elseif setprobe > 16
               pname = strrep(name,'.mat',sprintf('A.p%d.mat',ch));
            else
               pname = strrep(name,'.mat',sprintf('.p%d.mat',ch));
            end
%
%  Need to change this to load multiple file if necessary

            dvfile = [spkdir '/' strrep(probes(id(1)).filename,'t0.','dvdt.')];
            if length(id) > 1 & length(setprobe) == 1
                [a,sid] = sort([probes(id).first]);
                sid = id(sid);
                Spk.times = [];
                Spk.values = [];
                Spk.codes = [];

                for j = 1:length(sid)
                    filename = [dp '/Spikes/' probes(sid(j)).filename];
                    if exist(filename,'file')
                        a = load(filename);
                        chname = probes(sid(j)).var;
                        if isempty(Spk.times)
                            Spk = a.(chname);
                        else
                            Spk.times = [Spk.times; a.(chname).times];
                            Spk.codes = [Spk.codes; a.(chname).codes];
                            Spk.values = [Spk.values; a.(chname).values];
                            if size(a.(chname).times,1) > size(a.(chname).values,1)
                                fprintf('Some Missing Spike values in %s\n',filename);
                                if logfid
                                    fprintf(logfid,'Some Missing Spike values in %s\n',filename);
                                end
                            end
                        end
                    else
                        fprintf('No file %s\n',filename);
                    end
                end
                spkch = 'Spk';
            else
                a = load(pname);
                chid = (strmatch('Ch',fieldnames(a)));
                chnames = fieldnames(a);
                Chspk = a.(chnames{chid(1)});
                spkch = 'Chspk';
            end
                clear('a');

        elseif state.nospikes == 0
            if setprobe < 0 %forces rebuilding of probe list
                setprobe = -setprobe;
            end
            [a,b] = splitpath(name);
            if ignoreSpikeO == 2
                d = dir([b '/Expt*Times.mat']);
            else
                d = dir(b);
            end
                
            spkdir = [b '/Spikes'];
            if exist(spkdir,'dir') && ignoreSpikeO ~= 2
                d = dir(spkdir);
            else
                spkdir = b;
            end
            Chspk = [];
            for j = 1:length(d)
                if regexp(d(j).name,'.p[0-9][0-9,t]*.mat') %% d(j).bytes < 525000000 & np < 2
                    fprintf('loading %s\n',d(j).name)
                    a =  load([spkdir '/' d(j).name]);
                    chnames = fieldnames(a);
                    chnames = chnames(strmatch('Ch',chnames));
                    if ~isempty(chnames) %is data in file
                        np = np+1;
                        if strncmp(a.(chnames{1}).title,'4Trode',6)
                            probes(np).probe = sscanf(a.(chnames{1}).title,'4Trode%d');
                        else
                            probes(np).probe = sscanf(a.(chnames{1}).title,'Spike %d');
                        end
                        probes(np).var = chnames{1};
                        probes(np).filename = d(j).name;
                        probes(np).pathname = spkdir;
                        probes(np).first = a.(chnames{1}).times(1);
                        probes(np).last = a.(chnames{1}).times(end);
                        probes(np).nspk = length(a.(chnames{1}).times);
                        probes(np).traces = a.(chnames{1}).traces;
                        if strfind(d(j).name,'A.p')
                        probes(np).source = 2;
                        else
                        probes(np).source = 1;
                        end
                        if probes(np).probe == setprobe
                            Chspk = a.(chnames{1});
                        end
                        if length(a.(chnames{1}).times) > size(a.(chnames{1}).values,1)
                            nerr=nerr+1;
                            errs{nerr} = [d(j).name 'Missing some values'];
                            if state.showerrs
                            msgbox(errs{nerr},'APlaySpkFile Error!!','modal');
                            end
                            fprintf('%s\n',errs{nerr});
                        end
                    end
                elseif regexp(d(j).name,'.p[0-9][0-9]*.mat') %% d(j).bytes < 525000000 & np < 2
                    fprintf('loading %s\n',d(j).name)
                    a =  load([spkdir '/' d(j).name]);
                    chnames = fieldnames(a);
                    chnames = chnames(strmatch('Ch',chnames));
                    if ~isempty(chnames) %is data in file
                        np = np+1;
                        probes(np).probe = sscanf(a.(chnames{1}).title,'Spike %d');
                        probes(np).var = chnames{1};
                        probes(np).filename = d(j).name;
                        probes(np).first = a.(chnames{1}).times(1);
                        probes(np).last = a.(chnames{1}).times(end);
                        if probes(np).probe == setprobe
                            Chspk = a.(chnames{1});
                        end
                        if length(a.(chnames{1}).times) > size(a.(chnames{1}).values,1)
                            nerr=nerr+1;
                            errs{nerr} = [d(j).name 'Missing some values'];
                            msgbox(errs{nerr},'APlaySpkFile Error!!','modal');
                            fprintf('%s\n',errs{nerr});
                        end
                    end
                elseif regexp(d(j).name,'NewClusterTimes.mat') %% Ignore these
                elseif regexp(d(j).name,'ClusterTimes.mat') %% d(j).bytes < 525000000 & np < 2
                    suffix = 0;
                    csuffix = 0;
                    id = regexp(name,'\.[0-9,a]*\.mat');
                    if length(id)
                        suffix = sscanf(name(id(1)+1:end),'%d');
                        id = regexp(d(j).name,'Expt.[0-9,a]*');
                        if length(id)
                            csuffix= sscanf(d(j).name(id(1)+4:end),'%d');
                        end
                    end
%shouldn't need to load autoclustertimes except special circumstances
%since its put into ClusterTimes. When would you combine without making
%ClusterTimes? 
% - when the autoclusters are good!  Forget this can happen. 
%Only load autoclustertimes if clusertimes doesn't exist
                    if ~isempty(strfind((d(j).name),'AutoCluster'))
                        s = strrep(d(j).name,'AutoCluster','Cluster');
                        if sum(strcmp(s,{d.name}))
                            go = 0;
                        else
                            go = 1;
                        end
                    else
                        go = 1;
                    end
                    if suffix == csuffix && go
                    fprintf('loading %s\n',d(j).name);
                    clusterdate = d(j).datenum;
                    cname = [spkdir '/' d(j).name];
                    te = now;
                    a =  load(cname);
                    a.filename = cname;
                    if quickload == 0
                        dname = strrep(cname,'ClusterTimes','ClusterTimesDetails');
                        if exist(dname)
                            b = load(dname);
                            gotempty = 0;
                            for k = 1:length(b.ClusterDetails)
                                if isfield(b.ClusterDetails{k},'t')
                                    a.Clusters{k}.times = b.ClusterDetails{k}.t;
                                    a.Clusters{k}.xy = b.ClusterDetails{k}.xy;
                                    a.Clusters{k}.clst = b.ClusterDetails{k}.clst;
                                    %                            a.Clusters{k}.Evec = b.ClusterDetails{k}.Evec;
                                else
                                    gotempty = gotempty+1;
                                end
                            end
                        end
                        if gotempty || length(b.ClusterDetails) < length(a.Clusters) || gotempty
                            dname = strrep(cname,'ClusterTimes','AutoClusterTimesDetails');
                            if exist(dname)
                                b = load(dname);
                                for k = 1:length(a.Clusters)
                                    if isfield(b.ClusterDetails{k},'t') && ~isfield(a.Clusters,'xy')
                                        a.Clusters{k}.times = b.ClusterDetails{k}.t;
                                        a.Clusters{k}.xy = b.ClusterDetails{k}.xy;
                                        a.Clusters{k}.clst = b.ClusterDetails{k}.clst;
                                        %                            a.Clusters{k}.Evec = b.ClusterDetails{k}.Evec;
                                    end
                                end
                            end
                        end
                    else
                    end
                    nspkt = nspkt+1;
                    SpkTimes(nspkt).loadtime = (now-te)*60 *60*24;
                    id = strfind(d(j).name,'Expt');
                    if id
                        SpkTimes(nspkt).expno = sscanf(d(j).name(id(1)+4:end),'%d');
                    else
                        SpkTimes(nspkt).expno = 0;
                    end
                    SpkTimes(nspkt).section = 0;
                    id = regexp(name,'\.[0-9]*a\.mat');
                    if length(id)
                        SpkTimes(nspkt).section = 1;
                    end
                    
                    if strfind((d(j).name),'AutoCluster')
                        AutoClusters = a.Clusters;
                        SpkTimes(nspkt).auto = 1;
                    else
                        SpkTimes(nspkt).auto = 0;
                    end
                    SpkTimes(nspkt).datenum = clusterdate;
                    SpkTimes(nspkt).filename = d(j).name;
                    mint = 10e12;
                    maxt = 0;
                    for k = 1:length(a.Clusters)
                        if ~isempty(a.Clusters{k})
                            if a.Clusters{k}.auto == 0 && SpkTimes(nspkt).auto == 0
                                a.Clusters{k}.auto = 0;
                            else
                                a.Clusters{k}.auto = 1;
                            end
                            mint = min([mint min(a.Clusters{k}.times)]);
                            maxt = max([mint max(a.Clusters{k}.times)]);
                        end
                        if quickload == 0 && ~isfield(a.Clusters{k},'xy') && isfield(SpkTimes(1),'Clusters') &&isfield(SpkTimes(1).Clusters{k},'xy')
                            a.Clusters{k}.xy = SpkTimes(1).Clusters{k}.xy;
                            a.Clusters{k}.times = SpkTimes(1).Clusters{k}.times;
                            a.Clusters{k}.clst = SpkTimes(1).Clusters{k}.clst;
                        end
                    end
                    SpkTimes(nspkt).Clusters = a.Clusters;
                    if isfield(a,'FullVData')
                        SpkTimes(nspkt).FullVData = a.FullVData;
                    elseif isfield(b,'FullVData')
                        SpkTimes(nspkt).FullVData = b.FullVData;
                    end
                    SpkTimes(nspkt).trange = [mint maxt];
                    end
                    
                elseif regexp(d(j).name,'ClusterTimesDetails.mat') %% d(j).bytes < 525000000 & np < 2
                end
            end
            if isempty(Chspk) & exist('chnames','var') && ~isempty(chnames)
                Chspk = a.(chnames{1});
            end
            if ~isempty(probes)
                [a, id] = sort([probes.probe]);
                probes = probes(id);
                ip = unique([probes.probe]);
%sort by time of first spike keep track of total # spikes
                for j = ip;
                    id = find([probes.probe] == j);
                    [t, sid] = sort([probes(id).first]);
                    probes(id) = probes(id(sid));
                    nspk = 1;
                    for k = id;
                        probes(k).firsti = nspk;
                        nspk = nspk + probes(k).nspk;
                    end
                end
                save(idxfile,'probes');
            end
            clear a;
            spkch = 'Chspk';
                else %state.nospikes > 0
        [a,b] = splitpath(name);
        d = dir([b '/*ClusterTimes.mat']);
        end

    end
    if ~isempty(probes)
        [a, id] = sort([probes.probe]);
        probes = probes(id);
    end
    vnames = {'Ch30' 'Ch31' stimch};
    vlabels = {'Text' 'Events' 'Stim ON/OFF'};
    if state.needframes
        vlabels = [vlabels 'Frames'];
        vnames = [vnames framechname];
    end
    if state.nospikes == 0
        vlabels = [vlabels 'Spikes'];
        vnames = [vnames spkch];
    end
    for j = 1:length(vnames)
        missing(j) = ~exist(vnames{j},'var');
    end
    if sum(missing)
        msgbox(sprintf('%s Missing %s',name,vlabels{find(missing)}),'APlaySpkFile Error!!','modal');
        fprintf('%s Missing %s\n',name,vlabels{find(missing)});
        if logfid
            fprintf(logfid, '%s Missing %s\n',name,vlabels{find(missing)});
        
        fclose(logfid);
        end
        return;
    end
    Text = Ch30;
    if state.nospikes
        Spks.loaded = 0;
        Spks.times = [];
        Spks.codes = [];
        Spks.probe = 0;
    elseif ignoreSpikeO == NaN && length(probes) > 1   && Oprobe > 1
        for j = 1:length(probes)
            spkch = probes(j).var;
            p = probes(j).probe;
            AllSpikes{p} = eval(['CleanSpikes(' spkch ');']);
            AllSpikes{p}.times = AllSpikes{p}.times + timeoffset./10000;
        end
        Spks = AllSpikes;
    elseif ~eval(['isempty(' spkch ')']);
        if savedvdt & dvfile
            Spks = eval(['CleanSpikes(' spkch ',''dvfile'',dvfile);']);
        else
            Spks = eval(['CleanSpikes(' spkch ');']);
        end
        Spks.times = Spks.times+timeoffset./10000;
        clear(spkch);
    else
        if nspkt
            expts = unique([SpkTimes.expno]);
            expts = expts(expts > 0);
            np = 0;
            for j = expts
                id = find([SpkTimes.expno] == j & [SpkTimes.auto] == 0);
                aid = find([SpkTimes.expno] == j & [SpkTimes.auto] == 1);
                if length(aid) & length(id)

                    Clusters = SpkTimes(id(1)).Clusters;
                    for c = 1:length(SpkTimes(aid(1)).Clusters);
                        if c > length(Clusters) || isempty(Clusters{c})
                            Clusters{c} = SpkTimes(aid(1)).Clusters{c};
                        end
                        if ~isfield(Clusters{c},'excludetrialids')
                            Clusters{c}.excludetrialids{1} = [];
                        else
                            c = c;
                        end
                    end
                    if length(aid) == 2 && length(id) ==2 %Expt was split up
                        bClusters = SpkTimes(id(2)).Clusters;
                        for c = 1:length(SpkTimes(aid(2)).Clusters);
                            if c > length(bClusters) || isempty(bClusters{c})
                                Clusters{c}.times  = cat(2,Clusters{c}.times, SpkTimes(aid(2)).Clusters{c}.times);
                                Clusters{c}.clst  = cat(1,Clusters{c}.clst, SpkTimes(aid(2)).Clusters{c}.clst);
                                Clusters{c}.xy  = cat(1,Clusters{c}.xy, SpkTimes(aid(2)).Clusters{c}.xy);
                            else
                                Clusters{c}.times  = cat(2,Clusters{c}.times, bClusters{c}.times);
                                Clusters{c}.clst  = cat(1,Clusters{c}.clst, bClusters{c}.clst);
                                Clusters{c}.xy  = cat(1,Clusters{c}.xy, bClusters{c}.xy);
                                if isfield(bClusters{c},'excludetrialids')
                                    Clusters{c}.excludetrialids{1} = cat(2,Clusters{c}.excludetrialids{1}, bClusters{c}.excludetrialids);
                                end
                                    
                            end
                        end
                    end
                elseif length(aid)
                    Clusters = SpkTimes(aid).Clusters;
                else
                    Clusters = SpkTimes(id).Clusters;
                end
                np = max([np length(Clusters)]);
                EClusters{j} = Clusters;
                suffixlist(j) = j;
            end
            for j = 1:np
                Spks{j}.times = [];
                Spks{j}.codes = [];
                Spks{j}.cx = [];
                Spks{j}.cy = [];
            end
            for k = 1:length(EClusters)
                Clusters = EClusters{k};
            for j = 1:length(Clusters)
                if Clusters{j}.shape == 0
                    Clusters{j}.sign = 1;
                    Clusters{j}.crit = 0;
                end
                if ~isfield(Clusters{j},'sign') || Clusters{j}.sign == 0
                    Clusters{j}.sign = 1;
                end
                if isfield(Clusters{j},'times')
                    probes(j).probe = j;
                    if isempty(Clusters{j}.times)
                        Spks{j}.times = [];
                    elseif diff(size(Clusters{j}.times) > 1)
                        Spks{j}.times = cat(1,Spks{j}.times,Clusters{j}.times');
                    else
                        Spks{j}.times = cat(1,Spks{j}.times,Clusters{j}.times);
                    end
                    if quickload %just loaded times of classified spikes
                        Spks{j}.codes(:,1) = ones(size(Clusters{j}.times));
                        for c = 1:length(Clusters{j}.next)
                            if isfield(Clusters{j}.next{c},'times')
                            x = [Spks{j}.codes(:,1)' ones(size(Clusters{j}.next{c}.times)) * (1+c)];
                            Spks{j}.codes = x';
                            Spks{j}.times = cat(1,Spks{j}.times,Clusters{j}.next{c}.times');
                            end
                        end
                    else
                        if isfield(Clusters{j},'clst')
                            Spks{j}.codes(:,1) = Clusters{j}.clst-1;
                        else
                            Spks{j}.codes(:,1) = Clusters{j}.xy(:,1).*Clusters{j}.sign > Clusters{j}.crit.*Clusters{j}.sign;
                        end
                        Spks{j}.cx = Clusters{j}.xy(:,1);
                        Spks{j}.cy = Clusters{j}.xy(:,2);
                    end
%                    Spks{j}.codes(:,2) = Spks{j}.codes(:,1);
%mahal 1 9s 2D, mahal 4 in 1D, mahal 3 id ND, 0 if not used
                    if isfield(Clusters{j},'fitdprime')
                        Spks{j}.dips = [Clusters{j}.fitdprime(1) Clusters{j}.mahal(4) Clusters{j}.mahal(1)];
                        if Clusters{j}.space(1) == 6 && Clusters{j}.mahal(3) > Clusters{j}.mahal(1)
                            Spks{j}.dips(3) = Clusters{j}.mahal(3);
                        end
                    elseif isfield(Clusters{j},'mahal')
                        Spks{j}.dips = [NaN Clusters{j}.mahal(4) Clusters{j}.mahal(1)];
                    elseif isfield(Clusters{j},'dipsize')
                        Spks{j}.dips = [Clusters{j}.hdip(1) Clusters{j}.dipsize(1) abs(Clusters{j}.dprime(1))];
                    else
                        Spks{j}.dips = [NaN NaN NaN];
                    end
                    if isfield(Clusters{j},'dropi')
                        Spks{j}.dropi = Clusters{j}.dropi(3);
                    end
                    Spks{j}.next = [];
                    if isfield(Clusters{j},'next') && iscell(Clusters{j}.next)
                        for k = 1:length(Clusters{j}.next)
                            if isfield(Clusters{j}.next{k},'fitdprime')
                                Spks{j}.next{k}.dips(1) = Clusters{j}.next{k}.fitdprime(1);
                            end
                            if isfield(Clusters{j}.next{k},'mahal')
                                Spks{j}.next{k}.dips(2) = Clusters{j}.next{k}.mahal(4);
                                Spks{j}.next{k}.dips(3) = Clusters{j}.next{k}.mahal(1);
                            end
                            if isfield(Clusters{j}.next{k},'dropi')
                                Spks{j}.next{k}.dropi = Clusters{j}.next{k}.dropi(3);
                            end
                        end
                    end
                    Spks{j}.mahal = [Clusters{j}.mahal(1) Clusters{j}.mahal(2)];
                    Spks{j}.suffix = k; 
                    Spks{j}.crit = Clusters{j}.crit;
                    Spks{j}.clustersign = Clusters{j}.sign;
                    if isfield(Clusters{j},'excludetrialids')
                        Spks{j}.excludetrialids{1} = Clusters{j}.excludetrialids;
                    else
                        Spks{j}.excludetrialids{1} = [];
                    end
                    if isfield(Clusters{j},'missingtrials')
                        Spks{j}.excludetrialids{1} = union(Spks{j}.excludetrialids{1},Clusters{j}.missingtrials);
                    end
                    if length(Spks{j}.excludetrialids{1}) > 0
                        a = max(Spks{j}.excludetrialids{1});
                    end
                end
            end
            end
        else
            Spks.times = 0;
            Spks.codes = [0 0 0 0];
            Spks.values = 0;
        end
    end
    Events = Ch31;
    idxfile = strrep(name,'.mat','idx.mat');
    if idxfile & ~exist(idxfile,'file')
        mkidx = 1;
    elseif mkidx < 0
        mkidx = 0;
    elseif idxfile % idx file must exist - check its date
            d = dir(idxfile);
            dd = datenum(d.date);
            d = dir(name);
            if datenum(d.date) > dd % Expt file newer - rebuild
                fprintf('Rebuilding index: %s newer than %s\n',name,idxfile);
                mkidx = 1;
            end
    end
  
end


if isstruct(name)
    if isfield(name,'title') & isfield(name,'values') &... 
            size(name.values,2) == 46
        Spks = name;
        name = 'Unnamed'
    end
end

forcefix = 0;
argon = {};
j = 1;
while j <= nargin-1
    vg = varargin{j};
    if isstruct(vg) 
        if isfield(vg,'text')
            Text = vg;
        elseif isfield(vg,'codes');
            Events = vg;
        end
    elseif strncmpi(vg,'alltrials',8)
        argon = {argon{:} varargin{j}};
    elseif strncmpi(vg,'Defaults',4)
        j = j+1;
        defaults = varargin{j};
        if isfield(defaults,'starttrial')
            starttrial = defaults.starttrial;
        end
    elseif strncmpi(vg,'fixlfp',5)
        if strncmpi(vg,'fixlfpforce',8)
            forcefix = 1;
        end
        lfpfile = strrep(name,'.mat','A.lfp.mat');
        fixfile = strrep(name,'.mat','.lfp.mat');
% with 8 channels, everything is in one file - there is no 'A.lfp.mat' 
        if ~exist(lfpfile,'file') && exist(fixfile,'file')
            lfpfile = fixfile;
            forcefix = 1;
        end
        if exist(lfpfile,'file') && (~exist(fixfile,'file') || forcefix) ...
                && exist(mainsname,'var')
            load(lfpfile);
            if isfield(LFP.Header,'MainsNoise')
                fprintf('%s Already fixed\n',lfpfile);
                if logfid
                fprintf(logfid,'%s Already fixed\n',lfpfile);
                end
            else
            [LFP, avgs, NoiseAmp] = FixLFPMains(LFP,mainsch.times .* 10000);
            LFP.Header.amps = LFPGains(LFP);
            a = LFP.Header.amps ./ max(LFP.Header.amps);
            nch = sum(LFP.Header.chanlist > 0);
            if std(a(find(a > 0.1))) > 0.2 && nch <= 8  %% 8 channel probe tends to have mixed LFP gain
                LFP.Header.needscale = 1;
            else
                LFP.Header.needscale = 0;
            end
            save(fixfile,'LFP','NoiseAmp');
           if logfid > 0
               fprintf(logfid, '%s Fixed LFP in %s\n',datestr(now),fixfile);
           end
            end
        else
            fprintf('No LFP file %s\n',lfpfile);
           if logfid > 0
               fprintf(logfid, 'No LFP file %s\n',lfpfile);
           end
            
        end
        if strncmpi(vg,'fixlfponly',8)
            return;
        end
            
    elseif strncmpi(vg,'findtrial',5)
        j = j+1;
        argon = {argon{:} varargin{j-1} varargin{j}};
        findtrial = varargin{j};
    elseif strncmpi(vg,'mkufl',4)
        MkUfl(name, Ch30);
    elseif strncmpi(vg,'name',4)
        j = j+1;
        name = varargin{j};
   elseif strncmpi(vg,'noidx',5)
       mkidx = 0;
    elseif strncmpi(vg,'online',4)
        onlinedata = 1;
    elseif strncmpi(vg,'cluster',2)
        j = j+1;
        thecluster = varargin{j};
    elseif strncmpi(vg,'play',4)
        playspikes = 1;
    elseif strncmpi(vg,'rfs',3)
        MkUfl(name, Text,'overwrite');
    elseif strncmpi(vg,'relistonly',8)
        mkidx = 2;
    elseif strncmpi(vg,'relist',4)
        mkidx = 1;
    end
    j = j+1;
end

Header.Name = BuildName(name);
Header.loadname = name; 

Events.times = Events.times * 10000;
%spk times need to be ints for trigsdf.
if ~isempty(Spks)
    if iscell(Spks)
        for j = 1:length(Spks)
            if isfield(Spks{j},'times')
                Spks{j}.times = round(Spks{j}.times * 10000);
            end
        end
    else
        Spks.times = round(Spks.times * 10000);
    end
end
if isstruct(UstimV) && ~isempty(UstimV)
UstimV.times = round(UstimV.times * 10000);
end
Text.times = Text.times * 10000;
[rfstr, rfdat] = MkUfl(name, Text);

if strncmp(Text.comment,'GridData',8) || strncmp(Text.comment,'uProbe',8)
    Expt.DataType = Text.comment;
else
    Expt.DataType = 'Spike2';
end
Expt.Header.DataType = Expt.DataType;

if ~mkidx
    load(idxfile);
    if ~isfield(Expt,'errs')
        Expt.errs = {};
    end
    if ~isfield(Expt,'DataType')
        Expt.DataType = 'Spike2';
    end
% make sure name is is correct windows/unix form
    Expt.Header.Name = BuildName(Expt.Header.Name);
    Expt.Header.loadname = Header.loadname;
    Expt.Header.rfstr = rfstr;
    Expt.Header.rf = rfdat;
    Expt.setprobe = setprobe;
    Expt.Header.DataType = Expt.DataType;
    if isfield(Expt,'Comments') & isfield(Expt.Comments,'Peninfo');
        txt = GetPenInfo(Ch30);
        txt = strrep(txt,'Contact CNT','ContactCNT');
        Expt.Comments.Peninfo.trode = txt;
        id = strfind(txt,'Contact');
        if length(id)
            x = id(1);
            id = strfind(txt(x:end),' ');
            x = sscanf(txt(id(1)+x:end),'%d');
            Expt.Comments.Peninfo.probesep = x;
        end
    end
    if iscell(Spks)  && isfield(Spks{1},'values')
        AllData.AllSpikes = Spks;
        Spks = AllData.AllSpikes{1};
    elseif iscell(Spks)
        for j = 1:length(Spks)
            if isempty(Spks{j})
            else
                AllData.AllClusters(j) = Spks{j};
                for k = 1:nspkt
                    if isfield(SpkTimes,'FullVData') && ~isempty(SpkTimes(k).FullVData)
                        AllData.FullVData = SpkTimes(k).FullVData;
                    end
                end
            end
        end
        
        AllData.datenum = clusterdate;
        AllData.quickload = quickload;
    else
        AllData.Spikes = Spks;
    end
    if nspkt
        Expt.ClusterLoadTimes = [SpkTimes.loadtime];
    end
    Raw.stimlvl = stimlvl;
    Raw.stimch = stimch;
    if isstruct(UstimV)
        AllData.UstimV = UstimV;
    end
    if timeoffset
        Events.times = Events.times + timeoffset;
        f = {'Start' 'End' 'estimes' 'bstimes' 'TrueEnd' 'FalseStart'}
        for j = 1:length(f)
            if isfield(Expt.Trials,f{j})
                Expt.Trials.(f{j}) = Expt.Trials.(f{j})+timeoffset;
            end
        end
    end
       
    AllData.Events = Events;
    Expt.Probes = probes;
    if Expt.Header.Spike2Version < 1.23
        Expt = FixExpt(Expt,'ed');
    end
    iExpts = Expts;
    if isfield(Expt.Trials,'ve') && length(Expt.Trials.ve) < length(Expt.Trials.End)
        ve = mean(Expt.Trials.ve);
        Expt.Trials.ve(end:length(Expt.Trials.End)) = ve;
    end
    if isfield(Expts,'firsttrial')
            if timeoffset
        argon = {argon{:} 'timeoffset' timeoffset};
            end
            
            if ~isfield(Expts,'e3')
                [Expts.e3] = deal('e0');
            end
   [Expts, Expt] = SortExpts(Expts, Expt.Trials, Expt.Header, thecluster, Expt, state, argon{:});
   Expts = AddComments(Expts,Expt);
    if ~exist('Exptlist','var')
        ExptList = MkExList(Expts);
        newlist = 0;
    else
        newlist = 1;
    end
   if Expt.newerrs || newlist  %% why not separate ifs? why rewrite idxfile, when newerrs?
       WriteErrors(idxfile, Expt);
       tExpts = Expts;
       Expts = iExpts;
       if timeoffset == 0 %don't re-write the idxfile
       save(idxfile,'Expt','Expts','ExptList');
       end
       Expts = tExpts;
   end
    end
   if logfid >= 0
       fclose(logfid);
   end
   Expt.state = state;
   if saveexpts
       SaveExpts(name, Expts);
   end
    return;
end
Expt.state = state;
if 0 & idxfile & exist(idxfile,'file')
    load(idxfile);
    return;
end

    if iscell(Spks)  && isfield(Spks{1},'values')
        AllData.AllSpikes = Spks;
        Spks = AllData.AllSpikes{1};
    elseif iscell(Spks)
        for j = 1:length(Spks)
            if isempty(Spks{j})
            else
                AllData.AllClusters(j) = Spks{j};
            end
        end
        Spks = Spks{1};
    else
    AllData.Spikes = Spks;
    end



nt = 0;
nx = 0;
Expts = [];

%id = find(ismember(Spks.codes(:,1),thecluster));
Expt.setprobe = setprobe; 
frametimes = [];
bstimes = [];
estimes = [];

if exist(framechname,'var')
    if isfield(framech,'level')
        id = find(framech.level == 1);
        frametimes = framech.times(id) .* 10000
    else
    frametimes = eval([framechname '.times * 10000']);
    end
    Header.frameperiod = median(diff(frametimes));
else
    Header.frameperiod = 167;
end
if ~isempty(ustimmarkch)
    ustimmarkch.times = ustimmarkch.times .* 10000;
else
    ustimes = [];
end


if fixup(1) && stimlvl.level(end) == 1 && stimlvl.level(1) == 0
    stimlvl.times = stimlvl.times(1:end-1);
    stimlvl.level = stimlvl.level(1:end-1);
end

if exist(stimch,'var') & isfield(stimlvl,'level') & length(stimlvl.times) > 1
%    id = find(Ch8.level == 1);
% need to add an extra event at the end of each list to avoid issues with
% empty finds on the last trial. But make it long after to avoid any
% confusion with the real one;
    maxt = max([max(Text.times) max(Events.times) max(stimlvl.times)]);

    if isfield(stimlvl,'inverted') && stimlvl.inverted
        bstimes = stimlvl.times(stimlvl.level == 0) * 10000;
        bstimes(end+1) = maxt+50000;
        estimes = stimlvl.times(stimlvl.level == 1) * 10000;
        estimes(end+1) = maxt+50010;
    else
        bstimes = stimlvl.times(stimlvl.level == 1) * 10000;
        bstimes(end+1) = maxt+50000;
        estimes = stimlvl.times(stimlvl.level == 0) * 10000;
        estimes(end+1) = maxt+50010;
    end
end


bsid = strmatch('bs',Text.text);
bsstimes = Text.times(bsid);
esid = strmatch('es',Text.text);
esstimes = Text.times(esid);

tstore = zeros(size(bsid));
tstore(strmatch('bss',Text.text(bsid,:)))=1;
Trials.bsstimes = bsstimes;
Trials.esstimes = zeros(size(bsstimes));
for j = 1:length(esstimes)
    id = find(bsstimes < esstimes(j));
    Trials.esstimes(id(end)) = esstimes(j);
end
%Trials.esstimes(id) = NaN;

tic;
opid = strmatch('op',Text.text);
for j = 1:length(opid)
    str= sscanf(Text.text(opid(j),:),'op%d');
    if ~isempty(str)
        storing(j) = str;
    else
        storing(j) = 0;
    end
end
storing = bitand(storing,STOREBIT);
%storeonoff is a list of times where storing is toggled
storeonoff = 1+find(abs(diff(storing)) > 0); 
if isempty(storeonoff)
    storeonoff = [1 1];
end

ve = strmatch('BGCS Version',Text.text);
if ve
    version = sscanf(Text.text(ve(1),:),'BGCS Version %f');
else
    version = 1.1;
end

ids = strmatch('fz',Text.text);
if ids
fzs = textscan(Text.text(ids,:)','fz%n');
framerate = mean(fzs{1});
%    framerate = sscanf(Text.text(ids(1),:),'fz%n');
else
    if isfield(defaults,'fz')
        framerate = defaults.fz;
    else
        framerate = 96;
    end
end

ids = strmatch('nf',Text.text);
if ids
    fstr = Text.text(ids,1:5);
    fstr(fstr==0) = ' ';
    fzs = textscan(fstr','nf%n');
    nomdur = prctile(fzs{1},90) .* 10000/framerate;
end
instim = 0;
inexpt = 0;
nextonoff = 1;
if onlinedata
    storestate = 1;
else
    storestate = storing(1);
    tonoff = Text.times(opid(storeonoff(nextonoff)));
end

if onlinedata
Events.store = ones(size(Events.times));
else
tic;
Events.store = zeros(size(Events.times));
if ~storing(storeonoff(1)) & storing(1) %% first event is an off
    onid = find(Events.times < Text.times(opid(storeonoff(1))));
    Events.store(onid) = 1;
end

for j = 1:length(storeonoff)
    if storing(storeonoff(j))
        if length(storeonoff) > j
        onid = find(Events.times >= Text.times(opid(storeonoff(j))) ...
            & Events.times < Text.times(opid(storeonoff(j+1))));
        else
        onid = find(Events.times >= Text.times(opid(storeonoff(j))));
        end
        Events.store(onid) = 1;
    end
end
fprintf('Store Index takes %.2f\n',toc);
end

if onlinedata && sum(Events.codes(:,1)==STARTEXPT) == 0
    fprintf('Missing Expt Start First code at %.1f  ',Events.times(1)./10000);
    if logfid > 0
        fprintf(logfid,'Missing Expt Start First code at %.1f  ',Events.times(1)./10000);
    end
    id = find(Text.codes(:,1) == 1 & Text.times > 0.1);
    if ~isempty(id)
        Events.times = [Text.times(id(1)); Events.times];
        Events.store = [1; Events.store];
        Events.codes = [[STARTEXPT 0 0 0 ]; Events.codes];
        fprintf('True Start at %.1f\n',Text.times(id(1))./10000);
    end
end
nonstore = 0;
lastend = 0;
ntrials = sum(Events.codes(:,1) == FRAMESIGNAL);
trynew = 0;
if trynew %didn't help
Trials.Start(1:ntrials) = 0;
Trials.End(1:ntrials) = 0;
Trials.Trial(1:ntrials) = 0;
Trials.TrueEnd(1:ntrials) = 0;
Trials.Startev(1:ntrials) = 0;
Trials.stored(1:ntrials) = 0;
Trials.Result(1:ntrials) = 0;
Trials.serdelay(1:ntrials) = 0;
Trials.bstimes(1:ntrials) = 0;
Trials.delay(1:ntrials) = 0;
Trials.endelay(1:ntrials) = 0;
Trials.estimes(1:ntrials) = 0;
Trials.id(1:ntrials) = 0;
Trials.FalseStart(1:ntrials) = 0;
end

%remove stimON going to off at start of trial. 
if estimes(1) < bstimes(1) && length(estimes) == length(bstimes)+1
    estimes = estimes(2:end);
elseif length(bstimes) == length(estimes)+1 %crash can leave trailing start e.g. ruf1914
    bstimes = bstimes(1:end-1);
end


isi = bstimes(2:end)-estimes(1:end-1);

id = find(isi < 400); %shouldnt happen, but see lemM209.5
if length(id)
    nerr = nerr+1;
    errs{nerr} = sprintf('Bads ISIs (%d), starting at %.0f)\n',length(id),bstimes(id(1)));
    if state.showerrs
        msgbox(errs{nerr},'APlaySpkFile Error!!','modal');
    end
    Result.badisi = id;
    fprintf('%d impossible isis\n',length(id));
% if this is a rapid pulse sequence, sometimes sent by binoc.
% then need to remove the first bstime preceding the first short isi
% N.B. this only works for first train
    if bstimes(id(1)+1)-bstimes(id(1)) < 200
        id = cat(1,id(1)-1,id);
    end
    useid = setdiff(1:length(bstimes), id+1);
    bstimes = bstimes(useid);
    estimes = estimes(useid);
    isi = bstimes(2:end)-estimes(1:end-1);
else
    Result.badisi = [];
end


readmethod = 0;
settrials = 0;
if state.method == 1
    badbad = [];
    if onlinedata
        id = strmatch('bs',Text.text);
    else
        id = strmatch('bss',Text.text);
    end
bid = id;
bsstimes = Text.times(id);
    if onlinedata
        id = strmatch('es',Text.text);
    else
        id = strmatch('ess',Text.text);
    end
esstimes = Text.times(id);
exendid = find((Events.codes(:,1) == ENDEXPT | Events.codes(:,1) == CANCELEXPT) & Events.store > 0);
exstartid = find(Events.codes(:,1) == STARTEXPT & Events.store > 0);

if length(bsstimes) < length(bstimes)
    t = bsstimes(1)-bstimes;
    id = find(t > -1000 & t  <10000);  %only 1 event from -500ms to + 100ms) from text mark
    fprintf('%d StimOns before first saved trial\n',id(end));
    bstimes = bstimes(id:end);
    estimes = estimes(id:end);
end
evid = [];
bsid = [];
allfsid = [];
ExptStart = Events.times(exstartid);
exendid = exendid(exendid > exstartid(1));
ExptEnd = Events.times(exendid);
ExptCode = Events.codes(exendid,1);
if isempty(ExptEnd)
    ExptEnd = Events.times(end);
    ExptCode = 2; %Crahsed out files treated as good
end
ExptStart(length(exstartid)+1) = ExptEnd(end) +1;

for j = 1:length(exstartid)
    ts = ExptStart(j);
    id = find(ExptEnd > ExptStart(j)  & ExptEnd < ExptStart(j+1));
    if isempty(id) %force end Expt if notthing - probabaly a crash
        if j < length(exstartid)
            te = ExptStart(j+1)-10000;
            Expts(j).result = ENDEXPT;
        else
        te = Events.times(end);
        Expts(j).result = ENDEXPT;
        end
    else
        te = ExptEnd(id(1));
        Expts(j).result = ExptCode(id(1));
    end
    id = find(Events.times > ts & Events.times < te);
    evid = cat(1,evid,id);
%take care of any hanging estime after expt end. Can happen with
%fixed/crashed files
    eid = find(estimes > ts & estimes < te);
    if isempty(eid)  %no trials in this expt
        Expts(j).result = -1;
    else
%if final trial in expt is a badfix, and its short, then final FRAMESIGNAL
%code in serial line can be AFTER the StimOn Channel goes to off.  So
%search for all FRAMESIGNALS before end Expt (no more esttimes anyway
%becuase of previous test.
    fsid = find(Events.codes(id,1) ==  5 & Events.times(id) < te);
    bid = find(bstimes > ts & bstimes < te);
    bids{j} = bid;
    fsids{j} = id(fsid);
    bsid = cat(1, bsid, bid);
    allfsid = cat(1, allfsid, id(fsid));
    Result.startcounts = [length(fsid) length(bid)];
    Result.name = name;
    if sum(Result.startcounts) == 0
        fprintf('No Trials for block %d\n',j);
        Result.bsdelay = [];
    elseif length(fsid) == length(bid)
        fprintf('Ex %d Stimlevel and FrameSignal Lengths match %d)\n',j,length(fsid));
        bsidx(id(fsid)) = bid;
        Result.bsdelay = Events.times(id(fsid)) - bstimes(bid);
    else
        fprintf('Length Mismatch for Stimlevel (%d) and FrameSignal (%d)\n',length(bid),length(fsid));
        FindMissingTimes(Events, Text, bstimes, estimes);
        Result.bsdelay = [];
    end
    end
    Expts(j).start = ts;
    Expts(j).end = te;
    if length(bid) && Expts(j).result >= 0
        Expts(j).starti = bid(1);
        Expts(j).endi = bid(end);
    else
        Expts(j).starti = NaN;
        Expts(j).endi = NaN;
    end
end
%Expts = Expts([Expts.result] ==ENDEXPT);
evid = unique(evid);
bsid = unique(bsid);
fsid = find(Events.codes(evid,1) ==  5); %frame signal
fstimes = Events.times(evid(fsid));
fsid = allfsid;
fstimes = Events.times(fsid);
for j = 1:length(Expts)
    Expts(j).firsttrial = find(bsid == Expts(j).starti);
    Expts(j).lasttrial = find(bsid == Expts(j).endi);
    good(j) =  ~isempty((Expts(j).firsttrial));
end

Expts = Expts(good);

if isempty(Expts)
    fprintf('No Expts in %s\n',name);
    return;
end
bfid = find(Text.codes(:,1) == BADFIX & Text.codes(:,4) ==2);
if length(fsid) == length(bsid)
   Result.bsdelay = fstimes - bstimes(bsid);
   readmethod = 1;
%Trials.EndTxt idetnifies the text event that matched the end of trial n   
  for j = 1:length(fsid)
      Trials.Start(j) = bstimes(bsid(j));
      Trials.End(j) = estimes(bsid(j));
      Trials.bstimes(j) = Trials.Start(j);
      Trials.estimes(j) = Trials.End(j);
      Trials.Trial(j) = j+starttrial;
      Trials.Result(j) = 1;
      id = find(esstimes > fstimes(j));
      if length(id) == 0
          Trials.Result(j) = -1;
     elseif j <= length(fsid)
         Trials.EndTxt(j) = esstimes(id(1));
         if j < length(fsid) && esstimes(id(1)) > fstimes(j+1)
             a = find(Text.times(bfid) < fstimes(j+1));
             if ~isempty(a) && Text.times(bfid(a(end))) > fstimes(j) %Bad fixation trial
                 Trials.Result(j) = 0;
                 Trials.EndTxt(j) = Text.times(bfid(a(end)));
             else
                 Trials.Result(j) = -1;
             end
         end
      else
         Trials.EndTxt(j) = 0;
      end
      if ~isempty(frametimes) % have VTR channel
                id = find(frametimes > Trials.Start(j) & frametimes < Trials.Start(j)+500);
                if ~isempty(id)
                    Trials.delay(j) = frametimes(id(1)) - Trials.Start(j);
                    Trials.Start(j) = frametimes(id(1));
                    Trials.FalseStart(j) = 0;
                else
                    Trials.FalseStart(j) = 1;
                    Trials.delay(j) = NaN;
                end

          
          id = find(frametimes > Trials.End(j));
          if ~isempty(id) & frametimes(id(1))-Trials.End(j) < 500
              Trials.endelay(j) = Trials.End(j) - frametimes(id(1));
              Trials.End(j) = frametimes(id(1));
              Trials.TrueEnd(j) = frametimes(id(1));
          else
          end
      end
  end
  %badfix is always (?? what about Sa?) detected by spike2, so the stimoff
  %marker will always be after this.
  bid = find(Events.codes(:,1) == 11); %bad fix
  for j = 1:length(bid)
      id = find(Trials.bstimes < Events.times(bid(j)));
      if length(id)
      Trials.Result(id(end)) = 0;
      bfdelay(id(end)) = Trials.estimes(id(end))-Events.times(bid(j));
      else
          badbad = [badbad Events.times(bid(j))];
      end
  end
  bid = find(Text.codes(:,1) == 11 & Text.codes(:,4) == 2); %bad fix from spike2
  gbid = strmatch('BAD Saccade',Text.text(bid,:));  %Did complete trial, but invalid psych sacc
  bid = setdiff(bid, bid(gbid));
  for j = 1:length(bid)
      id = find(Trials.bstimes < Text.times(bid(j)));
      if length(id)
      Trials.Result(id(end)) = 0;
      bfdelay(id(end)) = Trials.estimes(id(end))-Text.times(bid(j));
      else
          badbad = [badbad Text.times(bid(j))];
      end
  end
      


  settrials = 1;
  nt = length(fsid);
else
    badcount = [];
    for j =1:length(fsids)
        if length(fsids{j}) ~= length(bids{j})
            badcount(j) = j;
        end
    end
    fprintf('Cant use New Read method - mismatched counts at%s\n',sprintf(' %d',badcount))
    readmethod = -1;
end
end



if settrials == 0 %were not set by new method
    fprintf('Using old Read Method\n');
    readmethod = 0;
    id = find(Text.codes(:,1) == BADFIX);
    bftimes = Text.times(id);
    for j = 1:size(Events.codes,1)

% if storage turned of mid-stim, don't want to miss ENSTIM marker
    if Events.codes(j,1) ~= ENDSTIM
        storestate = Events.store(j);
    end        
    if Events.codes(j,1) == FRAMESIGNAL
        nowt = Events.times(j);
        if storestate
            if nt & Trials.Result(nt) < 0
                nt = nt;
            end
            nt = nt+1;
        Trials.Start(nt) = Events.times(j);
        Trials.End(nt) = Trials.Start(nt)+nomdur; %% just in case an online file is missing end
        if findtrial  & Trials.Start(nt) > findtrial
            findtrial = 0;
        end
        Trials.Startev(nt) = Events.times(j);
        Trials.stored(nt) = storestate;
        Trials.Result(nt) = 1;
        if nt > 1 & length(Trials.End) < nt-1
            nerr = nerr+1;
            errs{nerr} = sprintf('Missing end Trial %d (EX %.0f, start %.0f)\n',nt-1,inexpt,Trials.Start(nt-1));
            fprintf('Missing end Trial %d (EX %.0f, start %.0f)\n',nt-1,inexpt,Trials.Start(nt-1));
%            Trials.End(nt-1) = NaN;
        end
        Trials.Trial(nt) = nt+starttrial;
        instim = 1;
%the event time can be just before the stim signal since it is not delayed
%to the sync (? correct explanation. Definitely what happens
% is the times of teh StimulusON digital event markers
%frametimes are the Vsync digital markers (every frame)
% the digital step can be >20ms after the serial signal is received, e.g.
% in ruf1989 at 209.43 sec
% in lem017 at 120.6 it is nearl 80ms late. at 289.6 its 300ms late. Looks
% like we could check for immediately preceding STARTSTIM (6) being before 
% StimON to check for this

        id = find(bstimes < Events.times(j)+300);
        
        if Events.times(j) > 238170000
            bstimes(id(end));
        end
        if id
            Trials.serdelay(nt) = Events.times(j) - bstimes(id(end));
            Trials.bstimes(nt) = bstimes(id(end));
        end
        if ~isempty(id) & bstimes(id(end)) > lastend 
            if id(end) < length(bstimes) & bstimes(id(end)+1) - Trials.Start(nt) < 10 %< 1ms to next = probably early
                Trials.Start(nt) = bstimes(id(end)+1);
            elseif bstimes(id(end)) > lastend
                Trials.Start(nt) = bstimes(id(end));
            end
            if Trials.serdelay(nt) > 10000
                nt = nt;
            end
            Trials.stored(nt) = storestate;
%here, Trials.Start is the time of the Digital event marker
%If vertical retrace was recorded, set start time to next one of these
%To a first approximation, it is the find(framemtimes ......  calls (here
%and for ENDSTIM) that take all the time in this loop. 
            if ~isempty(frametimes) % have VTR channel
                id = find(frametimes > Trials.Start(nt) & frametimes < Trials.Start(nt)+500);
                if ~isempty(id)
                    Trials.delay(nt) = frametimes(id(1)) - Trials.Start(nt);
                    Trials.Start(nt) = frametimes(id(1));
                    Trials.FalseStart(nt) = 0;
                else
                    Trials.FalseStart(nt) = 1;
                    Trials.delay(nt) = NaN;
                end
            else
%                Trials.Start(nt) = Events.times(j);
                Trials.FalseStart(nt) = 2;
                Trials.delay(nt) = NaN;
            end
        else
            Trials.Start(nt) = Events.times(j);
            Trials.delay(nt) = NaN;
            if isempty(id)
                Trials.FalseStart(nt) = 1;
            elseif id(end) < length(bstimes) & bstimes(id(end)+1) - Events.times(j) < 800 ...
%                    & Events.codes(j-1) == STARTSTIM ... %< 1ms to next = probably early
            nerr = nerr+1;
                errs{nerr} = sprintf('StimON at %.2f is %.1f ms late but STARTSTIM at %.2f',...
                bstimes(id(end)+1),(bstimes(id(end)+1)-Events.times(j))./10,Events.times(j-1));
            fprintf('%s\n',errs{nerr}); 
                Trials.FalseStart(nt) = 0;
                
            else
%Serial input can be very late if Spike2 got busy. Use the DIO stimon -
%this is the true start
                Trials.FalseStart(nt) = Events.times(j) - bstimes(id(end));
                Trials.Start(nt) = bstimes(id(end));
                 nerr = nerr+1; errs{nerr} = sprintf('Missing StimON at %.2f %.2f), but STARTSTIM at %.2f',Events.times(j),bstimes(id(end)),Events.times(j-1));
                fprintf('%s\n',errs{nerr});
                id = find(frametimes > Trials.Start(nt));
                if ~isempty(id) 
                    Trials.delay(nt) = frametimes(id(1)) - Trials.Start(nt);
                    Trials.Start(nt) = frametimes(id(1));
                end
            end
            Trials.stored(nt) = storestate;
        end 
        else
            nonstore = nonstore+1;
        end %if storestate
  
    elseif Events.codes(j,1) == ENDSTIM & storestate & nt
    if abs(Events.times(j) - 62740381) < 200
            Trials.End(nt) = Events.times(j);
            Trials.endelay(nt) = NaN;
    end
        if estimes
            id = find(estimes < Events.times(j)+500);
            if(id)
                Trials.TrueEnd(nt) = estimes(id(end));
                Trials.End(nt) = estimes(id(end));
                Trials.endelay(nt) = NaN;
                Trials.estimes(nt) = estimes(id(end));
  %if this is out by 400ms, probably failed to find correct end mark
                if Trials.TrueEnd(nt) < Events.times(j) - 4000  
                    fprintf('End event %.3f but marker %.3f\n',...
                        Events.times(j)./10000,Trials.TrueEnd(nt)./10000);
                    if Trials.End(nt) < Trials.Start(nt) & length(estimes) > id(end) & ...
                        estimes(id(end)+1) - Events.times(j) < 10000
                        Trials.TrueEnd(nt) = estimes(id(end));
                        Trials.End(nt) = estimes(id(end));
                        Trials.endelay(nt) = NaN;
                        Trials.estimes(nt) = estimes(id(end));
                    end
                end

                if ~isempty(frametimes) % have VTR channel
                    id = find(frametimes > Trials.End(nt));
                    if ~isempty(id) & frametimes(id(1))-Trials.TrueEnd(nt) < 500
                        Trials.endelay(nt) = Trials.End(nt) - frametimes(id(1));
                        Trials.End(nt) = frametimes(id(1));
                        Trials.TrueEnd(nt) = frametimes(id(1));
                    else
                    end
                end
            else
                Trials.TrueEnd(nt) = 0;
            end
            
        end
        Trials.End(nt) = Events.times(j);
        Trials.Result(nt) = 1;
        if (Trials.End(nt) - Trials.Start(nt)) < 1000
            instim = 0;
        end
        instim = 0;
        if Trials.TrueEnd(nt)
            lastend = Trials.TrueEnd(nt);
        else
            lastend = Trials.End(nt);
        end
    elseif Events.codes(j,1) == ENDTRIAL & storestate
        if instim
%  can't figure this out here because the BADFIX is only recorded in text,
%  not SampleKey (becuase this is send from Spike2, not received by, and
%  setting codes for sample keys is such a pain. But maybe should make all
%  of these events with code2 set to indicate it is from Spike2?
% Seems like this happens when fixation is broken just BEFORE stimulus on,
% but Spike2 has not registered this yet e.g. ruf2000 at 8867.9
%So, check badfix times from text (bftimes). If there is a badfix before
%this event, but after Stim start, nothing to worry about.
           id = find(bftimes < Events.times(j) & bftimes > Trials.Start(nt));
           Trials.Result(nt) = -1;  % this will be set to 0 if a BadFix is found.
           if isempty(id) 
               bsid = find(bstimes < Events.times(j));
               esid = find(estimes < Events.times(j));
               if length(bsid) && length(esid)
                   dur = estimes(esid(end))-bstimes(bsid(end));
               else
                   dur = (Events.times(j) - Trials.Start(nt));
               end
               if dur > 1000
                   nerr = nerr+1;
                   errs{nerr} = sprintf('End Trial without End stim: %d (%.2f) dur %.1fms',...
                       nt-1,Events.times(j)/10000,(Events.times(j) - Trials.Start(nt))/10,dur/10);
                   fprintf('%s\n',errs{nerr});
                    Trials.Result(nt) = -2;  % this will be set to 0 if a BadFix is found.
               end
           end
            Trials.End(nt) = Events.times(j);
            Trials.TrueEnd(nt) = NaN;
        end
    elseif Events.codes(j,1) == BADFIX & storestate %% Doesn't happen. Badfix is in Text, because it is sent, not received
        Trials.End(nt) = Events.times(j);
        Trials.Result(nt) = 0;
        instim = 0;
    elseif Events.codes(j,1) == STARTEXPT & storestate
        if inexpt %close an existing expt (e.g. if crashed out)
            Expts(nx).end = Events.times(j);
            Expts(nx).lasttrial = nt;
        end
        nx = nx+1;
        Expts(nx).start = Events.times(j);
        Expts(nx).firsttrial = nt+1;
        inexpt = 1;
    elseif Events.codes(j,1) == STARTEXPT % non stored
        inexpt = 0;
    elseif Events.codes(j,1) == ENDEXPT & nx & storestate
        Expts(nx).end = Events.times(j);
        Expts(nx).lasttrial = nt;
        inexpt = 0;
        Expts(nx).result = ENDEXPT;
    elseif Events.codes(j,1) == CANCELEXPT & nx & storestate
        Expts(nx).end = Events.times(j);
        Expts(nx).lasttrial = nt;
        Expts(nx).result = CANCELEXPT;
        inexpt = 0;
    elseif Events.codes(j,1) == ENDEXPT
        nx = nx;        
    end
end

end
Trials.id = zeros(size(Trials.Trial))'; %needs to be a row 
if nt == 0
    Expts = [];
    return;
end
if inexpt
    Expts(nx).lasttrial =nt;
end
if length(Expts) > nx  && settrials == 0 %some Expts set in first pass never found
    fprintf('Only found %d of %d good Expts\n',nx,length(Expts));
end
ntrials = nt;
fprintf('Setting Trials Took %.2f\n',toc);
tic;

trial = 1;
k = 1;
ix = 1;
nx = 1;

%kludge. Should fill esstimes with when the badfix happened.
if readmethod == 1 && length(Trials.End) > length(Trials.esstimes)
    Trials.End = Trials.End(1:length(Trials.esstimes));
end
id = find(Trials.esstimes(1:length(Trials.End)) == 0);
Trials.esstimes(id) = Trials.End(id);

% trynew should work, but need to check with some xxx= strings.
%problem with new method was cell2str deblanked. mat2cell works better
%

trynew = 1;
if trynew
    tic;
    
aText.text = mat2cell(Text.text(:,1:end-1),ones(1,size(Text.text,1)),size(Text.text,2)-1);
id = strmatch('xxx=',Text.text);
ids = setdiff(1:size(Text.text,1),id);
aText.text = aText.text(ids);
xid = ones(size(id));
nl = size(Text.text,1)-length(id);
k=1;
m = 1;
for j = 1:size(Text.text,1)
    if strncmp(Text.text(j,:),'xxx=',4)
        k = k-1;
        aText.text{k} = [aText.text{k} Text.text(j,5:end-1)];
        xid(m) = j;
        m = m+1;
% if we ever go back to this, probably need deblank(Text.....
    end
    if length(aText.text{k}) == 0 
 % Blank lines get removed with the deblank that follows. This misaligns text and codes.
 % so macke sure lines aren't blank. can always find these lines later.
        aText.text{k} = 'blank';
    end
    k = k+1;
end

aText.times = Text.times(ids);
aText.codes = Text.codes(ids,:);
aText.text = deblank(aText.text);
fprintf('Reading Text took %.2f\n',toc);
end

for j = 1:length(Expts)
    Expts(j).midtrial = (Expts(j).lasttrial+Expts(j).firsttrial)/2;
end

ix = 1;
if trynew == 0
tic;
for j = 1:size(Text.text,1)
    if strncmp(Text.text(j,:),'xxx=',4)
        k = k-1;
        aText.text{k} = [aText.text{k} Text.text(j,5:end-1)];
    else
% if we ever go back to this, probably need deblank(Text.....
    aText.text{k} = Text.text(j,1:end-1);
    aText.times(k) = Text.times(j);
    aText.codes(k,:) = Text.codes(j,:);
    end
    if j == 206991
        k
    end
    if length(aText.text{k}) %%? safe may remove blank strings that have codes.
        k = k+1;
    else
 % Blank lines get removed with the deblank that follows. This misaligns text and codes.
 % so macke sure lines aren't blank. can always find these lines later.
        aText.text{k} = 'blank';
        k = k+1;
    end
end
aText.text = deblank(aText.text);
fprintf('Reading Text %.2f\n',toc);
end
%
%AddTxtFile allows problems with data files to be fixed by adding
%additional lines written by hand. Fornmat is
%time  text
%where time is an int in timestamp units (0.1ms)
AddTxtFile = strrep(name,'.mat','Add.txt');
fid = fopen(AddTxtFile,'r');
if fid > 0
    a = textscan(fid,'%d %s','delimiter','\n');
    fclose(fid);
    t = a{1};
    s = a{2};
    for j = 1:length(s)
        %?why do we look for whitespace?? removed Aug 2010.
        id = findstr(s{j},' ');
        id = [];
        if length(id)
            txt = s{j}(id(1)+1:end);
        else
            txt = s{j};
        end
        if t(j) < 0 %special case for fixing lines
            if strncmp(s{j},'cm=rf',5)
                id = strmatch('cm=rf',aText.text);
                for k = 1:length(id)
                    aText.text(id,:) = s(j);
                end
            end
        elseif t(j) == 0 || t(j) <= aText.times(1)
            aText.text = {txt aText.text{:}};
            aText.times = [0; aText.times];
            aText.codes = [0 0 0 0; aText.codes];
        else
            id = find(aText.times < t(j));
            id = id(end);
            aText.text = {aText.text{1:id} txt aText.text{id+1:end}};
            aText.times = [aText.times(1:id); t(j); aText.times(id+1:end)];
            aText.codes = [aText.codes(1:id,:); 0 0 0 0; aText.codes(id+1:end,:)];
        end
        if ~isempty(txt)
            Text.text(end+1,1:length(txt)) = txt;
        end
    end
end
Trials.Stimseq = {};
intrial = 0;
Peninfo.trode = '';
nrw = 0;
tic;
tstart = now;
lasttook = 0;
Stimulus.CorLoop = 0;
Stimulus.SpikeGain = 50; %default
Stimulus.id = 0;
gotend = 0;
txtid = [];
lastfix.fx = 0;
lastfix.fy = 0;
Stimulus.Flag = '';
fprintf('Text->stims .....',toc);
inexpt = 0;
Stimulus.OptionCode = '+se';
nfpj=0;
bsctr = 0;
for j = 1:length(aText.text)
    aText.text{j} =  deblank(aText.text{j});
    txt = aText.text{j};
    t = aText.times(j);

%    id = find(txt == 0);
%    txt = txt(1:id(1)-1);
    if length(txt)>1;
        ss = txt(1:2);
        if length(txt) > 2 & txt(3) == '='
            val= txt(4:end);
        else
            val= txt(3:end);
        end
        if aText.codes(j,4) == 2  % this was FROM spike 3
            if aText.codes(j,1) == 3 && instim == 1 %end stim
                instim = 2;
            end

        elseif aText.codes(j,1) == 5 %stim start
            instim = 1;
            Stimulus.Seedseq = {};  %% these must be set for each stim
            Stimulus.Stimseq = {};
            Stimulus.xoseq = {};
            Stimulus.yoseq = {};
            Stimulus.Phaseseq = [];
            Stimulus.cLseq = [];
            Stimulus.cRseq = [];
            gotend = 0;
            nfpj=0;
            bsctr = bsctr+1;
        elseif aText.codes(j,1) == 3 %end stim
            instim = 2;
        elseif aText.codes(j,1) == STARTEXPT %end stim
            inexpt = 1;
            if readmethod == 1
                id = find([Expts.end] < aText.times(j));
                ix = length(id)+1;
            end
        elseif aText.codes(j,1) == ENDEXPT %end stim
            if inexpt
                ix = ix+1;
            end
            inexpt = 0;
        elseif strncmp(txt,'{}',2) %bug!!
        txt = aText.text{j};
            
        elseif strncmp(txt,'EndStim',7) %finished reading all text related to last stim
            gotend = 1;
        elseif strncmp(txt,'exvals',6)
            a = sscanf(txt,'exvals %f %f %f %d');
            if isfield(Stimulus,'et')
                Stimulus.(Stimulus.et) = a(1);
            end
            if isfield(Stimulus,'e2')
                if strmatch(Stimulus.e2,'backMov')
                    Stimulus.(Stimulus.e2) = a(2);
                end
            end
            Stimulus.ex3val = a(3);
            if isfield(Stimulus,'e3')
                if strcmp(Stimulus.e3,'mixac')
                    Stimulus.(Stimulus.e3) = a(3);
                end
            end
        elseif strncmp(txt,'mixac',5)
            Stimulus.mixac = sscanf(txt(6:end),'%f');
        elseif strncmp(txt,'Off at',6) %Storage turned off - should be outside trial
            if intrial
                fprintf('Storage Off in Trial at %.1f',atext.times(j));
            end
            
        elseif strncmp(txt,'RightHemi',9) || strncmp(txt,'Electrode',8)
            txtid = [txtid j];
            Peninfo.trode = txt;
            a = InterpretLine(txt);
            Peninfo = CopyFields(Peninfo, a);
            id = strfind(txt,'Contact');
            if length(id)
                x = id(1);
                id = strfind(txt(id:end),' ');
                sscanf(txt(id+x:end),'%d',x);
                Peninfo.probesep = x;
            end

        elseif strncmp(txt,'cm=rf',5)
            a = sscanf(txt,'cm=rf%f,%f:%fx%f,%fdeg pe%f %f,%f fx=%f,fy=%f');
            Stimulus.rf = a;
        elseif strncmp(txt,'StartDepth',10)
            Stimulus.StartDepth = str2num(txt(11:end));
        elseif strncmp(txt,'CLOOP',5)
            Stimulus.CorLoop = 1;
%        elseif strncmp(txt,'id',2)
%           Stimulus.id = sscanf(txt(3:end),'%d')
        elseif strncmp(txt,'bt',2)
             g = sscanf(txt,'bt%d spkgain %f');
             if length(g) > 1 & g(2) > 1
                 Stimulus.SpikeGain = g(2);
             end
             if length(g) > 0
                 ExptStartTime = g(1);
             end
        elseif strncmp(txt,'rw',2)
            nrw = nrw+1;
            [Trials.rws(nrw), ok] = sscanf(val,'%f');
            Trials.rwset(nrw) = t;
        elseif strncmp(txt,'st',2)
            Stimulus.st = strmatch(val, stimnames,'exact');
            Stimulus.st = Stimulus.st -1;
        elseif strncmp(txt,'mtop=op',7)
            Stimulus.OptionCode = txt(8:end);
        elseif strncmp(txt,'fl+',3)
            Stimulus.Flag = txt(3:end);
        elseif strncmp(txt,'mtrP=',5)
            Stimulus.Phaseseq = sscanf(txt(6:end),'%d');
            if length(Stimulus.Phaseseq) > length(Stimulus.Stimseq)
                Stimulus.Phaseseq = Stimulus.Phaseseq(1:length(Stimulus.Stimseq));
            end
            if trial > length(Trials.Start)
                Trials.Phaseseq{trial} = Stimulus.Phaseseq;
            elseif instim ~= 1  && trial > 1 && t < Trials.Start(trial)
                Trials.Phaseseq{trial-1} = Stimulus.Phaseseq;
            end
            if Stimulus.id == 6136 || Stimulus.id > 570
                trial;
            end
        elseif strncmp(txt,'mtFn=',5) && length(txt) > 50
            framets = sscanf(txt(6:end),'%f');
            id = find(diff(framets(1:end-1)) > 1.5);
            if diff(framets(end-1:end)) > 2.8
                id = cat(1,id ,length(framets));
            end
            if instim ~= 1 && trial > 1
                if trial > length(Trials.Start)  ||  t < Trials.Start(trial)
                    Trials.rptframes{trial-1} = id;
                    Trials.ndrop(trial-1) = length(id);
                    Trials.framet{trial-1} = framets;
                elseif instim == 2 && t > Trials.End(trial)
                    Trials.framet{trial} = framets;
                    Trials.ndrop(trial) = length(id);
                    Trials.rptframes{trial} = id;
                end
                
            end
        elseif strncmp(txt,'mtrX=',5)
            Stimulus.xoseq = sscanf(txt(6:end),'%d');
            if instim ~= 1 && trial > 1
                if trial > length(Trials.Start)  ||  t < Trials.Start(trial)
                    Trials.xoseq{trial-1} = Stimulus.xoseq;
                elseif instim == 2 && t > Trials.End(trial)
                    Trials.xoseq{trial} = Stimulus.xoseq;
                end
                
            end
        elseif strncmp(txt,'mtrY=',5)
            Stimulus.yoseq = sscanf(txt(6:end),'%d');
            if instim ~= 1 && trial > 1
                if trial > length(Trials.Start)  ||  t < Trials.Start(trial)
                    Trials.yoseq{trial-1} = Stimulus.yoseq;
                elseif instim == 2 && t > Trials.End(trial)
                    Trials.yoseq{trial} = Stimulus.yoseq;
                end
            end
        elseif strncmp(txt,'mtco=',5)
            Stimulus.Stimseq = sscanf(txt(6:end),'%d');
        elseif strncmp(txt,'mtcL=',5)
            Stimulus.cLseq = sscanf(txt(6:end),'%x');
            if instim ~= 1 && trial > 1 
                if trial > length(Trials.Start) || t < Trials.Start(trial)
                    Trials.cLseq{trial-1} = Stimulus.cLseq;
                elseif instim == 2 && t > Trials.End(trial)
                    Trials.cLseq{trial} = Stimulus.cLseq;
                end
            end
        elseif strncmp(txt,'mtcR=',5)
            Stimulus.cRseq = sscanf(txt(6:end),'%x');
            if sum(Stimulus.cRseq < 0) > 1
                Stimulus.cRseq = sscanf(txt(6:end),'%x');
            end
            if instim ~= 1 && trial > 1 
                if trial > length(Trials.Start) || t < Trials.Start(trial)
                    Trials.cRseq{trial-1} = Stimulus.cRseq;
                elseif instim == 2 && t > Trials.End(trial)
                    Trials.cRseq{trial} = Stimulus.cRseq;
                end
            end
        elseif strncmp(txt,'mtrS=',5)
            if aText.codes(j,1) == ENDSTIM
                istim = 2;
            end
            Stimulus.Stimseq = sscanf(txt(6:end),'%d');
            if Stimulus.id >= 5430
                Stimulus.Stimseq;
            end
            if instim ~= 1 && trial > 1 
                if trial > length(Trials.Start) || t < Trials.Start(trial)
                    Trials.Stimseq{trial-1} = Stimulus.Stimseq;
                elseif instim == 2 && t > Trials.End(trial)
                    Trials.Stimseq{trial} = Stimulus.Stimseq;
                end
            else
                instim;
            end
        elseif strncmp(txt,'mtse=',5)
            Stimulus.Seedseq = sscanf(txt(6:end),'%d');
            if instim ~= 1 && trial > 1 
                if trial > length(Trials.Start)  ||  t < Trials.Start(trial)
                    Trials.Seedseq{trial-1} = Stimulus.Seedseq;
                elseif instim == 2 && t > Trials.End(trial)
                    Trials.Seedseq{trial} = Stimulus.Seedseq;
                end
            else
                instim;
            end
            
        elseif strncmp(txt,'Nf',2) %comes after end TRIAL, not every stim
            if ~instim & trial > 1
                Trials.Nf(trial-1) = str2num(val);
            end
        elseif strncmp(txt,'mtet=',5)
            id = strfind(txt,'Fr');
            Stimulus.Fr = sscanf(txt(id+3:end),'%d');
        elseif strncmp(txt,'mtxo=',5)
            if isfield(Expts,'firsttrial') & isfield(Expts,'lasttrial')
                lastix = ix;
            ix = find([Expts.midtrial] > trial);
            if ix
                ix = ix(1);
                Expts(ix).xovals = sscanf(txt(6:end),'%f');            
            end
            end
        elseif strncmp(txt,'mtei=',5)
            if isfield(Expts,'firsttrial') & isfield(Expts,'lasttrial')
                lastix = ix;
            ix = find([Expts.midtrial] > trial);
            if ix
                ix = ix(1);
                Expts(ix).e1vals = sscanf(txt(6:end),'%f');            
            else
                nerr = nerr+1;
                errs{nerr} = sprintf('No Expt for mtei at trial %d',trial);
                fprintf('%s\n',errs{nerr});
                ix = lastix;
            end
            end
        elseif strncmp(txt,'mte3',5)
            if isfield(Expts,'firsttrial') & isfield(Expts,'lasttrial')
            lastix = ix;
            ix = find([Expts.firsttrial] < trial+2 & [Expts.lasttrial] > trial);
            if ix
                Expts(ix).e3vals = sscanf(txt(6:end),'%f');
            else
                ix = lastix;
            end
            end
        elseif strncmp(txt,'mte2=',5)
            ix = FindExptn(Expts, readmethod, bsctr, trial, ix);
            Expts(ix).e2vals = sscanf(txt(6:end),'%f');
        elseif strncmp(txt,'Off at',5)
        elseif strncmp(txt,'EndExpt',5)
            if inexpt
                ix = ix+1;
            end
            inexpt = 0;
        elseif strncmp(txt,'sonull',5)
        elseif strncmp(txt,'NewConnect',7)
        elseif strncmp(txt,'BGCS Version',7)
        elseif strncmp(txt,'testflag',7)
        elseif strncmp(txt,'rptframes ',3)
            Stimulus.rptframes = sscanf(txt(10:end),'%d');
        elseif strncmp(txt,'ijump',5)
            Stimulus.ijump = sscanf(txt(6:end),'%d');
        elseif strncmp(txt,'seof',4)
            Stimulus.seof = sscanf(txt(5:end),'%d');
        elseif strncmp(txt,'/local',6)
            Stimulus.imprefix = txt;
        elseif strncmp(txt,'imve ',3)
            [a,b] = sscanf(txt,'imve %f,%f %f');
            Stimulus.imver = a(1);
            Stimulus.imseed = a(2);
            if length(a) > 2 & a(3) < 1
                Stimulus.impx = a(3);
            end
        elseif strncmp(txt,'Sa:',3)
            if aText.codes(j,4) == 1
            a =  aText.codes(j,1);
            end
        elseif strncmp(txt,'op',2)
            a = sscanf(txt(3:end),'%f,%f');
            if ~isempty(a)
                Stimulus.op = a(1);
            end
            if length(a) > 1
                Stimulus.optionb = a(2);
            end
            if isempty(Stimulus.op)
                fprintf('Missing op stim %d\n',trial);
                Stimulus.op = 0;
            end
        elseif strncmp(txt,'backMov',7)
            Stimulus.backMov = sscanf(txt(8:end),'%f');
        elseif strncmp(txt,'annTyp',6)
            Stimulus.annTyp = sscanf(txt(7:end),'%f');
        elseif ss(1) == '#'
            comment = ss;
        elseif ~isstrprop(ss(1),'alphanum')
            fprintf('Non-Printing Name %s\n',txt);
        elseif strncmp(txt,'vs',2) || strncmp(txt,'sq',2)
            try
                [a, ok] = sscanf(val,'%f');
                if ~ok
                    Stimulus.(ss) = val;
                else
                    Stimulus.(ss) = a(1);
                    if length(a) > 2
                        Stimulus.FlipDir = a(3);
                    end
                end
            catch
                ok = 0;
                Stimulus.FlipDir = 0;
            end
        elseif strncmp(txt,'fp',2)
            a = sscanf(val,'%f');
            if length(a) > 1
            if instim == 1
                nfpj = nfpj+1;
                Stimulus.dfx(nfpj) = a(1);
                Stimulus.dfy(nfpj) = a(2);
            else
                Stimulus.fx = a(1);
                Stimulus.dfx = a(1);                
                Stimulus.fy = a(2);
                Stimulus.dfy = a(2);
                nfpj = 0;
            end
            end
        elseif strncmp(txt,'fx',2) || strncmp(txt,'fy',2)
%if really in a stimulus, make note of new fx but keep original also
%if instim ==2, don't change anything.
            a = sscanf(val,'%f');
            if instim == 1
                Stimulus.(['d' ss]) = a;
            elseif instim == 0
                Stimulus.(ss) = a;
                Stimulus.(['d' ss]) = a;
                lastfix.(ss) = a;
            elseif instim == 2 %Post stim, but trial counter not yet incremented. Store value
                Stimulus.(ss) = a;
                Stimulus.(['d' ss]) = a;
                lastfix.(ss) = a;
            end
        else
%            if strmatch(ss,{'0' '1' '2' '3' '4'})
            if regexp(ss,'^[0-9]')
                ss = ['x' ss];
            end
            try
            [Stimulus.(ss), ok] = sscanf(val,'%f');
            if ~ok
                Stimulus.(ss) = val;
            end
            catch
                ok = 0;
            end
            if strncmp(ss,'et',2)
                Stimulus.(ss) = val;
                ix = FindExptn(Expts, readmethod, bsctr, trial, ix);
                if(ix) Expts(ix).et = val; end
            end
            if strncmp(ss,'e2',2) & length(ix)==1
                ix = FindExptn(Expts, readmethod, bsctr, trial, ix);
                Expts(ix).e2 = val;
            end
            if strncmp(ss,'e3',2) & length(ix) ==1
                ix = FindExptn(Expts, readmethod, bsctr, trial, ix);
                Expts(ix).e3 = val;
            end
        end
    end
    if aText.codes(j,1) == FRAMESIGNAL
            newstart = 1;
    else
        newstart = 0;
    end
    if (aText.codes(j,4) == 2 || isempty(txt))  && exist('Stimulus','var') %empty text = code from spike2 -> binoc
        correctdir = 0;
        if isfield(Stimulus,'OptionCode') && ...
                (~isempty(strfind(Stimulus.OptionCode,'+2a')) || ~isempty(strfind(Stimulus.OptionCode,'+afc'))) ...
                & isfield(Stimulus,'vs')
            [a,b] = max(abs([Stimulus.vs(1) Stimulus.sq(1)]));
%
% historically negative respdir means +ve sacccade value
% for exactly oblique saccades, sign of vertical component does it. 
            if b == 1
                correctdir = -sign(Stimulus.vs(1));
            else
                correctdir = -sign(Stimulus.sq(1));
            end
            Stimulus.rwdir = correctdir;
%            Stimulus.FlipDir = 1;
        end
        if trial == 3057
            instim = instim;
        end
%Real ON/Off times are set from the events above. But need to know that
%text following WURTZOK applies to the next stimulus. So instim = 2 means
%text has been received ending trial, but not officially over yet. 
%trial gets incremented at end stim. So response applies to trial -1
%but check that the time is sensible. Can get one of these events when
%storage is off, resetting the last stored trials
       if trial > 1
           tdelay = aText.times(j) - Trials.End(trial-1);
       else
           tdelay  = 0;
       end
       if isfield(Trials,'FlipDir')
           flipdir = Trials.FlipDir(trial-1);
       else
           flipdir = 1;
       end
        if aText.codes(j,1) == WURTZOKW & trial > 1 & tdelay < 10000
            Trials.RespDir(trial-1) = -1 * correctdir.*flipdir;
            Trials.Result(trial-1) = 1;
            instim = 2;
  %          Trials.score(trial-1) = 0;
 %           Trials.scoretime(trial-1) = aText.times(j);
        elseif aText.codes(j,1) == WURTZOK & trial > 1  & tdelay < 10000
            Trials.RespDir(trial-1) = 1 * correctdir.* flipdir;
            Trials.Result(trial-1) = 1;
            instim = 2;
 %           Trials.score(trial-1) = 1;
 %           Trials.scoretime(trial-1) = aText.times(j);
        elseif aText.codes(j,1) == BADFIX
            if trial <= length(Trials.End)
            if aText.times(j) > Trials.End(trial) && settrials== 0
                Trials.End(trial) = aText.times(j);
            end
            if t < Trials.Start(trial) 
%                if trial > 1 && t < Trials.End(trial-1)
%                    Trials.Result(trial-1) = 0;
%                end
            else
                Trials.Result(trial) = 0;
            end
            end
            instim = 2;
        else %shouldb't happen
            trial = trial;
       end
    end
    if readmethod == 1
        endtimes = Trials.EndTxt;
    else
        endtimes = Trials.End;
    end
    if length(endtimes) >= trial && t >= endtimes(trial) & instim & trial <= length(Trials.Start) && newstart == 0
%need to read past the end of the last trial a litle way to get things like
%Stimseq which come afterwards
        if trial > length(Trials.Start)
            break;
        end
        took = (now-tstart) * 24 * 60 *60;
        if took - lasttook > 30
            fprintf('%.0fsec..',took);
            lasttook = took;
        end
        if isfield(Stimulus,'st')
        Trials = SetTrial(Stimulus, Trials, trial, ntrials);
        Stimulus.CorLoop = 0;
        Stimulus.uf = '';
        Stimulus.fx = lastfix.fx;
        Stimulus.fy= lastfix.fy;
        Stimulus.rptframes = [];
        Stimulus.endevent = aText.codes(j,1);
        if ~isempty(ustimmarkch)
            marks = find(ustimmarkch.times > Trials.Start(trial)-0.1 & ustimmarkch.times < Trials.End(trial)+0.1 ...
                & ustimmarkch.codes(:,1) == 1);
            if length(marks)
            Trials.uStimt{trial} = ustimmarkch.times(marks);
            elseif isfield(Trials,'optionb') && bitand(Trials.optionb(trial),64)
            marks = find(ustimmarkch.times > Trials.Start(trial)-0.1 & ustimmarkch.times < Trials.End(trial)+0.1);
            end
            if bitand(Trials.optionb(trial),64)
            marks = find(ustimmarkch.times > Trials.Start(trial)-0.1 & ustimmarkch.times < Trials.End(trial)+0.1);
            end
        end
%        AllStimuli(trial) = Stimulus; %fails when extra element is added to Stimulus
        trial = trial+1;
        end
        instim = 0;
    elseif trial > length(Trials.Start)
        if gotend
          break;
        end
% if we have passed stim Start time, text must refer to next stimulus. 
%But don't increment based on codes coming back from Spike
    elseif t > Trials.Start(trial) && instim == 0 && aText.codes(j,4) ~= 2
        instim = 1;
    end
end


if isfield(Trials,'ve') && iscellstr(Trials.ve)
    for j = 1:length(Trials.ve)
        if strncmp(Trials.ve{j},'binoclean',8)
            x = sscanf(Trials.ve{j}(11:end),'%f');
            ve(j) = 10+x(1);
            if length(x) > 1
                ve(j) = ve(j) + x(2)./100;
            end
        else
            ve(j) = 0;
        end
    end
    Trials.ve = ve;
    if length(Trials.ve) < length(Trials.Start)
        Trials.ve(length(Trials.Start)) = median(ve);
    end
end

if isempty(Expts(end).et) && length(Expts) > 1
    Expts(end).et = Expts(end-1).et;
    Expts(end).e2 = Expts(end-1).e2;
    Expts(end).e3 = Expts(end-1).e3;
end

if ~isfield(Expts,'result') || isempty(Expts(end).result)
    Expts(end).result = 0;
end
cmid = strmatch('cm=',aText.text);
rfid = strmatch('cm=rf',aText.text);
bkid = strmatch('cm=noback',aText.text);
cmid = setdiff(cmid,[rfid; bkid]);
cmid = union(cmid,txtid);
Expt.Comments.text = {aText.text{cmid}};
Expt.Comments.times = aText.times(cmid);
Expt.Comments.Peninfo = Peninfo;

if length(Trials.op) < length(Trials.Start) 
        Trials = SetTrial(Stimulus, Trials, length(Trials.Start),ntrials);
end 
if isfield(Trials,'RespDir') &  length(Trials.RespDir) < length(Trials.Start)%fill in final trial
        Trials = SetTrial(Stimulus, Trials, length(Trials.Start),ntrials);
end
idx = find(Trials.Result < 0);
if ~isempty(idx)
    fprintf('%d Trials missing End\n',length(idx));
end
fprintf('Text->stims %.2f\n',toc);
tic;
fn = fieldnames(Trials);
ntrials = length(Trials.Start);
cellids = {};
for j = 1:length(fn)
    if iscell(Trials.(fn{j}))
        cellids = {cellids{:} fn{j}};
        if length(Trials.(fn{j})) < ntrials
            Trials.(fn{j}){ntrials} = '';
        end
    end
end
for j = 1:length(Trials.Start)
    for k = 1:length(cellids)
        if isempty(Trials.(cellids{k}){j})
            Trials.(cellids{k}){j} = '';
        end
    end
end
fprintf('Clearing %d cells took %.2f\n',length(cellids),toc);
trial = 1;
instim = 0;
colors = mycolors;
lastspk = 1;
maxspk = 1;
tpause = 0;
if playspikes
    GetFigure('SpikeV');
end
postdur = 500;

tic;
for trial = 1:length(Trials.Start)
    if isnan(Trials.End(trial)) && trial > length(Trials.End)
        Trials.End(trial) = Trials.Start(trial) + nomdur;
    end
    spkids = find(Spks.times > Trials.Start(trial)-preperiod & Spks.times < Trials.End(trial)+postperiod);
    Trials.Spikes{trial} = (Spks.times(spkids) - Trials.Start(trial));
    if spkids
        spkid(trial,:) = [spkids(1) spkids(end)];
        Trials.Cluster{trial} = Spks.codes(spkids,1);
        frames = [];
 %this must be just for caliration of timing. Can't see need for this
 %if we recording real data
        if frametimes & testframe
            for k = spkids'
                id = find(frametimes < Spks.times(k));
                if id
                    frames = [frames frametimes(id(end))-Spks.times(k)];
                end
            end
        end
        if isempty(frames)
            Trials.Frames(trial) = NaN;
        else
            Trials.Frames(trial) = frames(1);
        end
    else
        spkid(trial,:) = [0 0];
        Trials.Frames(trial) = NaN;
        Trials.Cluster{trial} = [];
    end
    if playspikes
        for spk = spkids; 
            adc = Spks.values(spk,:);
            plot(adc,'color',colors{Spks.codes(spk,1)+1});
            drawnow;
            hold on;
            energy(spk) = sum(diff(adc).^2);
            svar(spk) = var(adc);
            vw(spk) = svar(spk)/energy(spk);
            if tpause
                pause(tpause);
            end
        end
        %        GetFigure('SpikeXY');
        subplot(2,1,1);
        plot(energy(lastspk:maxspk),vw(lastspk:maxspk),'.');
        hold on;
        %       GetFigure('SpikeV');
        subplot(2,1,2);
        hold off;
    end
end
fprintf('Spikes  take %.3f\n',toc);
if exist(dstimch,'var') & exist('stimchange','var')
    for j = 1:length(stimchange.times)
        et = stimchange.times(j) * 10000;
        id = find(Trials.Start < et);
        if length(id)
            t = id(end);
            if Trials.End(t) > et
                id = find(frametimes < et);
                if ~isempty(id)
                    diffs(j) = et-frametimes(id(end));
                    et = frametimes(id(end));
                end
                Trials.Events{t}{1,2} = et - Trials.Start(t);
                Trials.Events{t}{1,1} = 'ns';
            end
        end
    end
else
    Trials.Events{trial} = [];
end


tic;
trial = 1;
for spk = 1:size(Spks.times)
    t = Spks.times(spk);
    if t > (Trials.End(trial)+postdur) & instim
        Trials.Spikes{trial} = (Spks.times(lastspk:maxspk) - Trials.Start(trial));
        spkid(trial,:) = [lastspk, maxspk];%
%        Trials(trial).Frames = frames;
        Trials.Cluster{trial} = Spks.codes(lastspk:maxspk,1);
        trial = trial+1;
        instim = 0;
        lastspk = spk;
        elseif t > Trials.Start(trial)
            instim = 1;
        end
        if instim & playspikes
            adc = Spks.values(spk,:);
            plot(adc,'color',colors{Spks.codes(spk,1)+1});
            energy(spk) = sum(diff(adc).^2);
            svar(spk) = var(adc);
            vw(spk) = svar(spk)/energy(spk);
            hold on;
        end
        maxspk = spk;
        if trial > length(Trials.Start)
            break;
        end
end
fprintf('Spikes Took %.3f\n',toc);
fixfile = strrep(idxfile,'idx','fix'); 
if exist(fixfile)
    load(fixfile);
    f = fields(fixTrials);
    [tid, fid] = ismember(fixTrials.id,Trials.id);
    tid = find(tid);
    for j = 1:length(f)
        if ~strcmp(f{j},'id')
            if ~isfield(Trials,f{j})
                Trials.(f{j}) = ones(size(Trials.Start));
            end
            Trials.(f{j})(fid(tid)) = fixTrials.(f{j})(tid);
        end
    end
end
%Trials = Trials(2:end);
AllData.datenum = clusterdate;
AllData.quickload = quickload;
AllData.Events = Events;
AllData.Text = aText;
AllData.Spikes = Spks;
if exist('Ch32','var') && strcmp(Ch32.title,'DigMark')
    Expt.DigMark = Ch32; %So its saved
end


stored = find(bitand([Trials.op],STOREBIT));
% bitand STOREBIT (16) identifies if storage was on;
%Trials = Trials(stored);

%Expts #20 seems to be empty, but storage was on...
if isfield(Trials,'TrueEnd') & length(Trials.TrueEnd) < length(Trials.End)
    Trials.TrueEnd(length(Trials.End)) = NaN;
end
if isfield(probes,'var')
for j = 1:length(probes)
   probes(j).var = strrep(probes(j).var,'A','');
end
end
Expt.Trials = Trials;
Expt.Spkid = spkid;
Expt.Probes = probes;
Expt.bstimes= bstimes;
Expt.estimes= estimes;
Expt.ExptList = Expts;
if exist('mainsch','var') && isfield(mainsch,'times')
    Expt.mainstimes = mainsch.times * 10000;
end

Header.Name = BuildName(name);
Header.Spike2Version = version;
Header.unstored = nonstore;
Header.CreationDate = CreationDate(Text);
Header.ReadMethod = readmethod;
if isfield(Expt,'DataType')
    Header.DataType = Expt.DataType;
end
clear Text;
if isempty(Expts) %nothing in this file
    fclose(logfid);
    return;
end

if mkidx == 1
        fprintf('Saving index %s\n',idxfile);
        if logfid > 0
        fprintf(logfid, '%s Saving index %s\n',datestr(now),idxfile);
        end
        Expt.bstimes = bstimes;
        Expt.estimes = estimes;
        Expt.Header = Header;
        if nerr > 0
        Expt.errs = errs;
        else
            Expt.errs = {};
        end
    save(idxfile,'Expt','Expts');
    iExpts = Expts;
end



if isfield(Expts,'firsttrial')
    if timeoffset
        args = {args{:} 'timeoffset' timeoffset};
    end
[Expts, Expt] = SortExpts(Expts, Trials, Header, thecluster, Expt, state);
Expts = AddComments(Expts,Expt);
end

if mkidx == 1
    tExpts = Expts;
    ExptList = MkExList(Expts);
    Expts = iExpts;
    Expt.starttrial = starttrial;
    save(idxfile,'Expt','Expts','ExptList');
    Expts = tExpts;
end
if saveexpts
    SaveExpts(name, Expts);
end
if mkidx
%    save(idxfile,'Expt','Expts');
end
fclose(logfid);
WriteErrors(idxfile, Expt);

function SaveExpts(name, Expts)
%make separate files for each expt on disk so that can access in parallel
%when needed
for j = 1:length(Expts)
    outfile = strrep(name,'.mat',['Expt' num2str(j) '.mat']);
    Expt = Expts{j};
    save(outfile,'Expt');
end


function ix = FindExptn(Expts, readmethod, bsctr, trial, ix)
if readmethod == 1
    trial = bsctr;
end
    
if isfield(Expts,'firsttrial') & isfield(Expts,'lasttrial')
    lastix = ix;
    ix = find([Expts.firsttrial] < trial+2 & [Expts.lasttrial] > trial);
    if ix
        ix = ix(end);
    else
        ix = lastix;
    end
end

function ExptList = MkExList(Expts)
ExptList = [];
if ~iscell(Expts)
    return;
end
for j =1:length(Expts)
        ExptList(j).expname = Expts{j}.Header.expname;
        ExptList(j).start = Expts{j}.Header.Start;
        ExptList(j).end = Expts{j}.Header.End;
        ExptList(j).et = Expts{j}.Stimvals.et;
        ExptList(j).e2 = Expts{j}.Stimvals.e2;
        ExptList(j).e3 = Expts{j}.Stimvals.e3;
end

function WriteErrors(idxfile, Idx)
if isfield(Idx,'errs') & length(Idx.errs)
    ename = strrep(idxfile,'idx.mat','err.txt');
    fid = fopen(ename,'a');
    fprintf(fid,'%s\n',Idx.errs{:});
    fclose(fid);
end

function Expts = AddComments(Expts, Expt)

if isfield(Expt,'Comments')

for j = 1:length(Expt.Comments.times)
    id = find(Expt.Trials.Start < Expt.Comments.times(j));
    if isempty(id)
        Expt.Comments.id(j) = Expt.Trials.id(1);
    else
        Expt.Comments.id(j) = Expt.Trials.id(id(1));
    end
end


    for j = 1:length(Expts)
        cid = find(Expt.Comments.times > Expts{j}.Header.trange(1)-100000 & ...
       Expt.Comments.times < Expts{j}.Header.trange(2)+100000);
        if ~isempty(cid)
            bid = strmatch('cm=back=',Expt.Comments.text(cid));
            for k = 1:length(bid)
                %these lines are before start expt. If after end its for
                %next expt
                if Expt.Comments.times(cid(bid(k))) < Expts{j}.Header.trange(2)
                    Expts{j} = ParseExptComment(Expts{j}, Expt.Comments.text{cid(bid(k))});
                end
            end
        end
        Expts{j}.Comments.text = {Expt.Comments.text{cid}};
        Expts{j}.Comments.times = {Expt.Comments.times(cid)};
        cid = strmatch('cm=VisualArea',Expt.Comments.text);
        for k = 1:length(cid)
            Expts{j}.Header.Area{k} = Expt.Comments.text{cid(k)}(15:end);
        end
    end
end



function Trials = SetTrial(Stimulus, Trials, trial, ntrials)

fn = fieldnames(Stimulus);
if isfield(Stimulus,'Ro') && isfield(Stimulus,'dx')
    ca = cos(Stimulus.Ro * pi/180);
    sa = sin(Stimulus.Ro * pi/180);
% need to recheck sign conventions here in replay...    
    Stimulus.Op = Stimulus.yo .* ca - Stimulus.xo .* sa;
    Stimulus.Pp = Stimulus.yo .* sa + Stimulus.xo .* ca;
    Stimulus.dO = Stimulus.dx .* sa - Stimulus.dy .* ca;
    Stimulus.dP = Stimulus.dy .* sa + Stimulus.dx .* ca;
end
for k = 1:length(fn)
    F = fn{k};
    if ~isfield(Trials,F)
%        Trials.F(1:ntrials) = NaN; % pre-allocate memory
    end
    if strcmp(F,'St')
        Trials.St(trial) = Stimulus.(F);
    elseif strncmp(F,'trode',5)
    elseif ischar(Stimulus.(F)) 
        if ~isfield(Trials,F) || ischar(Trials.(F)) || iscell(Trials.(F))
        if length(Stimulus.(F)) > 0
            Trials.(F){trial} = Stimulus.(F);
        end
        else
            fn{k}
        end
    elseif sum(strcmp(F,{'Seedseq' 'Stimseq' 'Phaseseq' 'cLseq' 'cRseq' 'xoseq' 'yoseq' 'rptframes'}))    
           Trials.(fn{k}){trial} = Stimulus.(fn{k});
    else
        if isempty(Stimulus.(F))
            Trials.(F)(trial) = NaN;
        else
            Trials.(F)(trial,1:length(Stimulus.(F))) = Stimulus.(F);
        end
    end
end
if isempty(Stimulus.st) || isnan(Stimulus.st) || isnan(Trials.st(trial))
    Stimulus.st
end
if isfield(Trials,'rwset') & isfield(Trials,'rws')
id = find(Trials.rwset < Trials.Start(trial));
eid = find(Trials.rwset < Trials.End(trial));
if length(id)
Trials.rw(trial) = Trials.rws(id(end));
end
end
if isfield(Trials,'RespDir') & length(Trials.RespDir) < length(Trials.Start)
    Trials.RespDir(length(Trials.Start)) = 0;
%    Trials.FlipDir(length(Trials.Start)) = 1;
end

function [Expts, Idx] = SortExpts(AllExpts, AllTrials, Header, thecluster, Idx,state,  varargin)
SpkDefs;
timeoffset = 0;
%stimnames = {'None', 'Gabor', 'RDS', 'Grating', 'bar', 'circle', 'rectangle', 'test', 'square', 'probe', '2grating', 'Cylinder', 'twobar', 'rls', 'annulus', 'rdssine', 'nsines'};
Expts = [];
spikid = Idx.Spkid;
if ~isfield(Idx,'newerrs')
    Idx.newerrs = 0;
end
if isfield(Header,'frameperiod')
frameperiod = Header.frameperiod;
else
frameperiod = 167;
end
findtrial = 0;
if state.alltrials
    usebadtrials = 1;
else
usebadtrials = 0;
end
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'alltrials',5)
        usebadtrials = 1;
    elseif strncmpi(varargin{j},'findtrial',5)
        j = j+1;
        findtrial = varargin{j};
    elseif strncmpi(varargin{j},'timeoffset',8)
        j = j+1;
        timeoffset = varargin{j};
    end
    j = j+1;
end

if isfield(AllTrials,'ve') && length(AllTrials.ve) < length(AllTrials.Start)
        AllTrials.ve(length(AllTrials.Start)) = median(AllTrials.ve);
end


for j = 1:length(AllExpts)
    if isempty(AllExpts(j).result)
        AllExpts(j).result = -1;
    end
end

AllExpts = AllExpts([AllExpts.result] >= 0);
STIM_BAR = 4;
STIM_GRATING = 3;
for j = 1:length(stimnames)
    if strmatch(stimnames{j},'bar')
        STIM_BAR = j-1; %starts with 0
    end
end
fn = fieldnames(AllTrials);
% make a list of fileds that are NOT automatically set 
ids = strmatch('Spikes',fn);
ids = [ids strmatch('OptionCode',fn)];
ids = [ids strmatch('Start',fn,'exact')];
ids = [ids strmatch('End',fn,'exact')];
ids = [ids strmatch('EndTxt',fn,'exact')];
ids = [ids strmatch('Result',fn)];
ids = [ids strmatch('op',fn)']; % has op and optionb
ids = [ids strmatch('Stimseq',fn)];
ids = [ids strmatch('Seedseq',fn)];
ids = [ids strmatch('xoseq',fn)];
ids = [ids strmatch('yoseq',fn)];
ids = [ids strmatch('rptframes',fn)];
ids = [ids strmatch('framet',fn)];
ids = [ids strmatch('cLseq',fn)];
ids = [ids strmatch('cRseq',fn)];
ids = [ids strmatch('Phaseseq',fn)];
ids = [ids strmatch('Cluster',fn)];
ids = [ids strmatch('Frames',fn)];
ids = [ids strmatch('StartEv',fn)];
ids = [ids strmatch('Events',fn)];
ids = [ids strmatch('rws',fn,'exact')];
ids = [ids strmatch('endelay',fn)];
ids = [ids strmatch('rwset',fn)];
ids = [ids strmatch('estimes',fn)];
ids = [ids strmatch('uStimt',fn)];
ids = [ids strmatch('Flag',fn,'exact')];
ids = [ids strmatch('bsstimes',fn,'exact')];
ids = [ids strmatch('esstimes',fn,'exact')];
ids = [ids find(strcmp('ex3val',fn))];
%ids = [ids strmatch('imver',fn)];
%ids = [ids strmatch('imse',fn)];

fn = fn(setdiff([1:length(fn)],ids));
for nf = 1:length(fn)
    if iscell(AllTrials.(fn{nf}))
    else
        if length(AllTrials.(fn{nf})) < length(AllTrials.Start)
            fprintf('Forcing end values for %s\n',fn{nf});
            AllTrials.(fn{nf})(length(AllTrials.Start)) = 0;
        end
    end
end
%non-zero values of falsestart should indicate the time gap
if isfield(AllTrials,'TrueStart') & ~isfield(AllTrials,'FalseStart')
    AllTrials.FalseStart = AllTrials.TrueStart;
    AllTrials.FalseStart(find(AllTrials.TrueStart == 0)) = 1;
    AllTrials.FalseStart(find(AllTrials.TrueStart == 1)) = 0;
end
if ~isfield(AllTrials,'Nf') || length(AllTrials.Nf) < length(AllTrials.Start)
    AllTrials.Nf(length(AllTrials.Start)) = NaN;
end

if ~isfield(AllTrials,'Events') || length(AllTrials.Events) < length(AllTrials.Start)
    AllTrials.Events{length(AllTrials.Start)} = [];
end
%If there are a minority of trials with False Starts (usually missing the
%StimON channel, or a long delay to the matching FRAMESIGNAL on the serial line
%set the result to 2, so that these can be exluded if necessary.
% only do this if the result is already not 0 (otherwise Add in Badfix
% trials
if isfield(AllTrials,'FalseStart') && sum(AllTrials.FalseStart > 0) < length(AllTrials.FalseStart)/5
    AllTrials.Result(find(AllTrials.FalseStart > 0 & AllTrials.Result > 0)) = 2;
end
    %fn is now a list of fields that are set for each trial automatically.
nexpts = 1;
    phasevals  = [0:360];
    phasevals(2:4) = [pi pi/4 3*pi/4];
    lphasevals = [0 pi pi 0]; %%see binoc.c SetRandomPhase();
    rphasevals = [0 pi 0 pi];
    for nx = 1:length(AllExpts)
    
    clear Trials;
    frpt = 1;
    a = AllExpts(nx).firsttrial:AllExpts(nx).lasttrial;
    if isfield(AllTrials,'sM') && median(AllTrials.sM(a)) == 26
        usebadtrials = 1;
    end
    if usebadtrials
        igood = 1:length(a);
        fn = {fn{:} 'Result'};
    else
        igood = find(AllTrials.Result(a) > 0);
    end
    nt = length(igood);
    igood = a(igood);
    nu = 0; %number of ustim pulses
    if nt <= 3 && AllExpts(nx).result ~= CANCELEXPT && state.showerrs
        msgbox(sprintf('%s Expt %d only %d trials',Header.Name,nx,nt),'APlaySpkFile Error!!','modal');
    end
    if ~isfield(AllExpts,'result')
        AllExpts(1).result = 1;
    end
    if nt> 3 & igood(1) < length(AllTrials.Start) & ismember(AllExpts(nx).result,[2 0])
       spkids = [];
       needfields = {};
       if isempty(AllExpts(nx).e3)
           fprintf('Empty Type Expt %d, trials %d - %d\n',nx, AllExpts(nx).firsttrial, AllExpts(nx).lasttrial)
       elseif strmatch(AllExpts(nx).e3,'ar')
           needfields = {needfields{:} 'wi' 'hi'};
        end
       for j = 1:nt
            if findtrial & AllTrials.Start(igood(j)) > findtrial
                findtrial = 0;
            end
            [Trials(j).Start] = AllTrials.Start(igood(j));
            [Trials(j).TrialStart] = AllTrials.Start(igood(j));
            if isfield(AllTrials,'TrueEnd') & AllTrials.TrueEnd(igood(j)) > 0
                [Trials(j).End] = AllTrials.TrueEnd(igood(j));
            else
                [Trials(j).End] = AllTrials.End(igood(j));
            end
            Trials(j).dur = Trials(j).End(end)-Trials(j).Start(1);
            if isfield(AllTrials,'optionb')
            [Trials(j).uStim] = bitand(AllTrials.optionb(igood(j)),64);
            end
            [Trials(j).op] = AllTrials.op(igood(j));
            [Trials(j).Trial] = AllTrials.Trial(igood(j));
            [Trials(j).id] = AllTrials.id(igood(j));
            if ~isempty(AllTrials.Events{igood(j)})
                Trials(j).Events = AllTrials.Events{igood(j)};
            end
            if isfield(Trials,'uStimt')
                nu = nu + legnth(Trials(j).uStimt);
            end
%            spkids = [spkids spkid(igood(j),1):spkid(igood(j),2)];
        end
 % cant do this. for some reason max spkids) is > size(Spikes.times)
%        Expt.gui.spks = unique(spkids);
%  problem is that we don't rebuilds Spkid if set s probe != 1, so save
%  Spkid no good
 %      
        for nf = 1:length(fn)
            if iscell(AllTrials.(fn{nf}))
                if strcmp(fn{nf},'uStimt')
                else
                nv = unique({AllTrials.(fn{nf}){igood}});
                Expt.Stimvals.(fn{nf}) = AllTrials.(fn{nf}){a(1)};
                end
            else
                nv = unique([AllTrials.(fn{nf})(igood)]);
            end
            
            if isnumeric(nv) & sum(~isnan(nv))> 1
                for j = 1:nt
                    Trials(j).(fn{nf}) = AllTrials.(fn{nf})(igood(j));
                end
                if strcmp(fn{nf},'Nf')
                    if sum(nv) == max(nv) % only 1 value + 0;
                        Expt.Stimvals.(fn{nf}) = max(nv);
                    end
                    Trials = rmfield(Trials,'Nf');
                end
                Expt.Stimvals.(fn{nf}) = prctile([AllTrials.(fn{nf})(igood)],50);
            elseif isnumeric(nv)
                Expt.Stimvals.(fn{nf}) = nv;
            end
        end
        if isfield(AllTrials,'dfx') && size(AllTrials.dfx,2) > 1
            for j = 1:nt
                id = find(AllTrials.dfx(igood(j),:) ~= 0);
                Trials(j).dfx = AllTrials.dfx(igood(j),id);
                Trials(j).dfy = AllTrials.dfy(igood(j),id);
            end
            
        end
        if isfield(AllTrials,'rf');
            Expt.Stimvals.rf = median(AllTrials.rf(igood,:));
            Expt.Stimvals.st = mode(AllTrials.st(igood));
        end
    duration = mean([Trials.End] - [Trials.Start]);
    et = Expt.Stimvals.et;
    e2 = Expt.Stimvals.e2;
    if isfield(Expt.Stimvals,'Fr') && Expt.Stimvals.Fr > 0
        frpt = Expt.Stimvals.Fr;
    end
    if isfield(AllExpts(nx),'e1vals') & ~isempty(AllExpts(nx).e1vals)
        Expt.e1vals = AllExpts(nx).e1vals;
    end
    if isfield(AllExpts(nx),'xovals') & ~isempty(AllExpts(nx).xovals)
        Expt.xovals = AllExpts(nx).xovals;
    elseif isfield(Expt,'xovals')
        Expt = rmfield(Expt,'xovals');
    end
    if strmatch('ce',{et e2})
        Expt.Stimvals.ce = median(abs(AllTrials.ce(igood)));
    end
     if strmatch(et,{'Op' 'Pp'},'exact') & isfield(Expt.Stimvals,'rf')
            ca = cos(Expt.Stimvals.rf(5) * pi/180);
            sa = sin(Expt.Stimvals.rf(5) * pi/180);
            rOp = Expt.Stimvals.rf(2) .* ca - Expt.Stimvals.rf(1) .* sa;
            rPp = Expt.Stimvals.rf(2) .* sa + Expt.Stimvals.rf(1) .* ca;
            Expt.Stimvals.rOp = rOp;
            Expt.Stimvals.rPp = rPp;
    if isfield(AllExpts(nx),'e1vals')
% don't mess with the values for interleaved extras
            sid = find(AllExpts(nx).e1vals > -1000);
            if strcmp(et,'Op')
                Expt.e1vals(sid) = AllExpts(nx).e1vals(sid) +rOp;
            elseif strcmp(et,'Pp')
                Expt.e1vals(sid) = AllExpts(nx).e1vals(sid) +rPp;
            end
        end
    end
    if isfield(AllExpts(nx),'e2vals') && ~isempty(AllExpts(nx).e2vals)
        Expt.e2vals = AllExpts(nx).e2vals;
        if strcmp(Expt.Stimvals.e2,'ce') && max(Expt.e2vals) > 1
            bid = find(Expt.e2vals > 1);
            err = sprintf('Fixing ce values (%s->1)',sprintf('%.1f ',Expt.e2vals(bid)));
            Idx = AddError(err, Idx, 0);
            Expt.e2vals(bid) = 1;
        end

    end
                if isfield(Expt.Stimvals,'Fs')
        xovals = [-16:16] .* Expt.Stimvals.Fs;
        yovals = [-16:16] .* Expt.Stimvals.Fs;
        end
    allev = [];
    serrid = [];
    nframes = [];
    psychtrial = 0;
    seqtrial = 0;
    crtrial = 0; %count # with contrast reversal
    timesexpt = 0;
    nu = 0; %count trials with uStimt;
    for k = 1:nt
        sid = find(ismember(AllTrials.Cluster{igood(k)}, thecluster));
        Trials(k).Spikes = round(AllTrials.Spikes{(igood(k))}(sid));
        Trials(k).count = sum(find(Trials(k).Spikes > 500 & Trials(k).Spikes < duration+500));
        Trials(k).sz = AllTrials.wi(igood(k));
        Trials(k).OptionCode = AllTrials.OptionCode{igood(k)};
        if isfield(AllTrials,'Flag')
        Trials(k).Flag = AllTrials.Flag{igood(k)};
        if strfind(Trials(k).Flag,'+mm')
            Trials(k).flatsurf = 1;
        else
            Trials(k).flatsurf = 0;
        end
        end
        if isfield(AllTrials,'uStimt')
        Trials(k).uStimt = round(AllTrials.uStimt{(igood(k))});
        nu = nu + length(Trials(k).uStimt);
        end
        if ~isempty(strfind(Trials(k).OptionCode,'+2a')) || ~isempty(strfind(Trials(k).OptionCode,'+afc'))
            psychtrial = psychtrial+1;
        end
        if strfind(Trials(k).OptionCode,'+fS')
            seqtrial = seqtrial+1;
        end
        if strfind(Trials(k).OptionCode,'+cr')
            crtrial = crtrial+1;
        end
        if strfind(Trials(k).OptionCode,'+x2')
            timesexpt = timesexpt+1;
        end
        if ~isnan(AllTrials.Frames(igood(k)))
            Trials(k).Frames = AllTrials.Frames(igood(k));
        end
        if bitand(LMONOC,Trials(k).op)
            Trials(k).me  = -1;
        elseif bitand(RMONOC,Trials(k).op)
            Trials(k).me  = 1;
        else
            Trials(k).me  = 0;
        end
        if isfield(Trials,'Fr')
            frpt = Trials(k).Fr;
        end
% for image seqs, Stimseq records the order of seeds, and there will not
% be a conversion to stimulus type
        if isfield(AllTrials,'Seedseq') && length(AllTrials.Seedseq{igood(k)}) > 1
            Trials(k).Seedseq = AllTrials.Seedseq{igood(k)};
        end
        for f = 1:length(needfields)
            Trials(k).(needfields{f}) = AllTrials.(needfields{f})(igood(k));
        end
        if isfield(AllTrials,'Stimseq') && length(AllTrials.Stimseq{igood(k)}) > 1 && ...
                ((isfield(Expt,'e1vals') && ~isempty(Expt.e1vals)) || isfield(Expt,'xovals')) 
            if AllTrials.Stimseq{igood(k)}(end) > length(Expt.e1vals) %last value sometimes junk               
               AllTrials.Stimseq{igood(k)} = AllTrials.Stimseq{igood(k)}(1:end-1);
            end
            if length(Expt.e1vals) && ~isempty(AllTrials.Stimseq{igood(k)}) && ...
                ((max(AllTrials.Stimseq{igood(k)}) < length(Expt.e1vals) && seqtrial > k/2) || strcmp(et,'backMov'))
            evid = AllTrials.Stimseq{igood(k)}+1;
            if frpt > 1
                evid = evid(1:frpt:end);
            end
            Trials(k).st = ones(size(evid)) * Expt.Stimvals.st;
            Trials(k).ce = ones(size(evid)) * Expt.Stimvals.ce;
            Trials(k).me = ones(size(evid)) * Trials(k).me(1);
            if strmatch(et,'Dc')
                ev = Expt.e1vals(evid);
                if frpt > 1
                    ev = ev(1:frpt:end);
                end
                Trials(k).(et) = AllTrials.Dc(igood(k));
                Trials(k).(Expt.Stimvals.e2)= ev;
                Trials(k).(e2)(find(ev == ISIGNALFRAME)) = AllTrials.(e2)(igood(k)); %% blanks
                if strcmp(Expt.Stimvals.e2,'or')
                    Trials(k).ori = AllTrials.(e2)(igood(k));
                end
            elseif strmatch(et,'backMov')
                ev = AllTrials.Stimseq{igood(k)}+1;
                Trials(k).(et) = ev;
             %   eb = Expt.e2vals(AllTrials.Stimseq{igood(k)}+1);
              %  Trials(k).(Expt.Stimvals.e2)= eb;
            else
                ev = Expt.e1vals(evid);
                Trials(k).(et) = ev;
                 if isfield(Expt,'e2vals') && length(Expt.e2vals) > max(AllTrials.Stimseq{igood(k)})
                    eb = Expt.e2vals(evid);
                    Trials(k).(Expt.Stimvals.e2)= eb;
                 else
                     eb = zeros(size(ev));
                end
            end
            if isfield(AllTrials,'xoseq') & length(AllTrials.xoseq{igood(k)}) > 0
                Trials(k).xo = xovals(AllTrials.xoseq{igood(k)}+1);
            end
            if isfield(AllTrials,'yoseq') & length(AllTrials.yoseq{igood(k)}) > 0
                Trials(k).yo = yovals(AllTrials.yoseq{igood(k)}+1);
            end
            if isfield(AllTrials,'cLseq') & length(AllTrials.cLseq{igood(k)}) > 0
                Trials(k).cL = AllTrials.cLseq{igood(k)};
                Trials(k).cL(Trials(k).cL < 0) = 0;
                Trials(k).cL(Trials(k).cL > 512) = 0;
            end
            if isfield(AllTrials,'cRseq') & length(AllTrials.cRseq{igood(k)}) > 0
                Trials(k).cR = AllTrials.cRseq{igood(k)};
                Trials(k).cR(Trials(k).cR < 0) = 0;
                Trials(k).cR(Trials(k).cR > 512) = 0;
            end
            if isfield(AllTrials,'Phaseseq') & length(AllTrials.Phaseseq{igood(k)}) > 0
                if Expt.Stimvals.sM ==13 && Expt.Stimvals.st == STIM_BAR
                Trials(k).lph = lphasevals(AllTrials.Phaseseq{igood(k)}+1);
                Trials(k).rph = rphasevals(AllTrials.Phaseseq{igood(k)}+1);
                elseif Expt.Stimvals.sM ==13 && Expt.Stimvals.st == STIM_GRATING
                Trials(k).lph = lphasevals(AllTrials.Phaseseq{igood(k)}+1);
                Trials(k).rph = rphasevals(AllTrials.Phaseseq{igood(k)}+1);
                else
                Trials(k).ph = phasevals(AllTrials.Phaseseq{igood(k)}+1);
                end
            end
            Trials(k).Start = Trials(k).Start + [0:length(ev)-1]' .* frameperiod * frpt;
            Trials(k).End = Trials(k).Start + frameperiod * frpt;
            Trials(k).st(find(ev == IBLANK)) = 0; %% blanks
            Trials(k).ce(find(ev == IUNCORR)) = 0; %% uncorr
            Trials(k).me(find(ev == ILEFTMONOC)) = -1; 
            Trials(k).me(find(ev == IRIGHTMONOC)) = 1;
                if size(Trials(k).(et),2) > 1
                    size(evid)
                end
                if isfield(Trials,e2) & size(Trials(k).(e2),2) > 1
                    size(evid)
                end
%in the Lopos X Ropos expt, could exceed max stim combinatinos, then mtei
%and mte2 strings only gave the inidividual Lpos/Rpos values, not all
%combinations
        elseif isfield(AllTrials,'Stimseq') && AllTrials.ve(end) < 4.85 && ...
                ~isempty(AllTrials.Stimseq{igood(k)}) && ...
                ~isempty(strfind(Trials(k).OptionCode,'+fS')) && ...
             max(AllTrials.Stimseq{igood(k)}) > length(Expt.e1vals)
         extras = sum(Expt.e1vals < -999);
         evid = AllTrials.Stimseq{igood(k)}-extras;
            Trials(k).st = ones(size(evid)) * Expt.Stimvals.st;
            Trials(k).ce = ones(size(evid)) * Expt.Stimvals.ce;
            Trials(k).me = ones(size(evid)) * Trials(k).me(1);
         n2 = length(Expt.e2vals)-extras;
         n1 = length(Expt.e1vals)-extras;
         e1 = mod(evid,n2)+extras+1;
         id = find(evid <0);
         e1(id) = evid(id)+extras+1;
         eb = floor(evid./length(Expt.e1vals))+1+extras;
         eb(id) = evid(id)+extras+1;
         Trials(k).(et) = Expt.e1vals(e1);
         Trials(k).(et)(id) = 0;
         ev = Expt.e1vals(e1);
         if ~isempty(Expt.e2vals)
         Trials(k).(e2) = Expt.e2vals(eb);
         end
            Trials(k).Start = Trials(k).Start + [0:length(ev)-1]' .* frameperiod * frpt;
            Trials(k).End = Trials(k).Start + frameperiod * frpt;
            Trials(k).st(find(ev == IBLANK)) = 0; %% blanks
            if Trials(k).id == 6138
            Trials(k).st(find(ev == IBLANK)) = 0; %% blanks
            end
            if isfield(AllTrials,'Phaseseq') & length(AllTrials.Phaseseq{igood(k)}) > 0
                if Expt.Stimvals.sM ==13 && Expt.Stimvals.st == STIM_BAR
                Trials(k).lph = lphasevals(AllTrials.Phaseseq{igood(k)}+1);
                Trials(k).rph = rphasevals(AllTrials.Phaseseq{igood(k)}+1);
                elseif Expt.Stimvals.sM ==13 && Expt.Stimvals.st == STIM_GRATING
                Trials(k).lph = lphasevals(AllTrials.Phaseseq{igood(k)}+1);
                Trials(k).rph = rphasevals(AllTrials.Phaseseq{igood(k)}+1);
                else
                Trials(k).ph = phasevals(AllTrials.Phaseseq{igood(k)}+1);
                end
            end
        elseif isfield(AllTrials,'Stimseq') && length(AllTrials.Stimseq{igood(k)}) > 2 && seqtrial > k/2
            if strcmp(Expt.Stimvals.et,'serange')
                if isfield(AllTrials,'xoseq') & length(AllTrials.xoseq{igood(k)}) > 0
                    Trials(k).xo = xovals(AllTrials.xoseq{igood(k)}+1);
                end
                if isfield(AllTrials,'yoseq') & length(AllTrials.yoseq{igood(k)}) > 0
                    Trials(k).yo = yovals(AllTrials.yoseq{igood(k)}+1);
                end
            end
            serrid = [serrid k];
            end
        end
        if isfield(AllTrials,'rptframes') && ~isempty(AllTrials.rptframes{igood(k)})
            Trials(k).rptframes = AllTrials.rptframes{igood(k)};
        end
        
        durs(k) = Trials(k).End(end) - Trials(k).Start(1);
        nframes(k) = length(Trials(k).Start);
    end
    if isfield(Trials,'CorLoop')
        id = find([Trials.CorLoop] == 0);
        if length(id)
            Trials = Trials(id);
            nframes = nframes(id);
            durs = durs(id);
        end
    end
    if isfield(Trials,'flatsurf')
        fs = unique([Trials.flatsurf]);
        if length(fs) == 1
            Trials = rmfield(Trials,'flatsurf');
            Expt.Stimvals.flatsurf = fs;
        end
    end

    if mean(nframes > 2)
       id = find(nframes > 1);
       Trials = Trials(id);
       bid = find(nframes  ==1);
       if length(bid) > 2
           fprintf('Removing %d Trials becuase no RC sequence\n',length(bid));
       end
    end
    if nu == 0  && isfield(Trials,'uStimt')
        Trials = rmfield(Trials,'uStimt');
    end
    if isfield(Trials,'uStim') && sum([Trials.uStim]) == 0
            Trials = rmfield(Trials,'uStim');
    end


    id = [];
    if isfield(Trials,'delay')
        id = find(~isnan([Trials.delay]));
    end
    if length(id) > length(Trials)/2
        nid = find(isnan([Trials.delay]));
        Expt.Trials = Trials(id);
        durs = durs(id);
        if length(nid)
            if isfield(Trials,'FalseStart')
            fprintf('Delay Nan at %.2f(%.2f), id%.0f\n',Trials(nid(1)).Start(1),Trials(nid(1)).FalseStart, Trials(nid(1)).id);
            else
            fprintf('Delay Nan at %.2f, id%.0f\n',Trials(nid(1)).Start, Trials(nid(1)).id);
            end
        end
    else
        Expt.Trials = Trials;
    end
    Expt.Header = Header;
    Expt.Header.trange(1) = AllExpts(nx).start;
    if isfield(AllExpts,'end')  && ~isempty(AllExpts(nx).end)
% If reaal End of expt is long time after last trial, probabaly means user did
% EndExpt. Still want to record real time of endexpt
        if length(AllExpts) < nexpts || Expt.Trials(end).End(end)+10000 > AllExpts(nx).end
            AllExpts(nx).end = Expt.Trials(end).End(end)+10000;
        else
            AllExpts(nx).end = AllExpts(nx).end+5000;
        end
    else %%online, unfinished, edpt
        if nx > 1  
            fprintf('Expt %d No end. Adding.\n');
        end
        AllExpts(nx).end = Expt.Trials(end).End(end)+10000;
    end
    if isfield(AllExpts,'end')  && ~isempty(AllExpts(nx).end)
        Expt.Header.trange(2) = AllExpts(nx).end;
    else
        Expt.Header.trange(2) = 0;
    end



    if isfield(Idx,'DigMark')
        t = Expt.Header.trange./10000;
        id = find(Idx.DigMark.times > t(1) & Idx.DigMark.times < t(2));
        if length(id)
            Expt.DigMark.times = Idx.DigMark.times(id);
            Expt.DigMark.codes = Idx.DigMark.codes(id,1);
        end
    end
    if psychtrial > nt/2
        Expt.Header.psych = 1;
    else
        Expt.Header.psych = 0;
    end
    if timesexpt > nt/2
        Expt.Stimvals.x2 = 1;
    else
        Expt.Stimvals.x2 = 0;
    end
    if seqtrial > nt/2
        Expt.Header.rc = 1;
    else
        Expt.Header.rc = 0;
    end
    Expt.Header.Options = '';
    if crtrial > nt/2
        Expt.Header.Options = [Expt.Header.Options '+cr'];
    end

    
    [a, expname, b, stimname] = Expt2Name(Expt);
    if isempty(expname)
        Expt.Header.expname = 'None';
    else
        Expt.Header.expname = [stimname '.' expname];
    end
    Expt.Header.Start = AllExpts(nx).start+timeoffset;
    Expt.Header.End = AllExpts(nx).end+timeoffset;
 
    if isfield(Expt.Trials,'cL')
 %       vals = [Expt.Trials.cL]; %unsafe syntax = varying Fr
        if median(nframes) > 1 %% subspace map for cL
            Expt.Header.expname = [Expt.Header.expname 'C'];
        end
    end
    if Expt.Header.rc
        Expt.Header.expname = [Expt.Header.expname 'RC'];
    end
    if min(durs) < 0
        id = find(durs < 0);
        err = sprintf('Expt %d at %.2f(id%d) has (%d) negative durations',nexpts,Expt.Trials(id(1)).Start(1),Expt.Trials(id(1)).id,length(id));
        Idx = AddError(err, Idx, state.showerrs);
        id = find(durs > 0);
        Expt.Trials = Trials(id);
    end
    if ~isfield(Expt,'e1vals') & isfield(Expt.Trials,et)
        Expt.e1vals = unique(cat(1,Expt.Trials.(et)));
    end

%replace off screen positiosn with 0, but only for 
%regular tuning curves
    if sum(strncmp(Expt.Stimvals.et,{'xo' 'yo' 'Op' 'Pp'},2)) && Expt.Header.rc == 0
        exv = [Expt.Trials.(Expt.Stimvals.et)];
        if size(exv,1) == 1
        id = find(abs(exv) > 35); %off screen
        for j = 1:length(id)
            Expt.Trials(id(j)).xo = Expt.Stimvals.rf(1);
            Expt.Trials(id(j)).yo = Expt.Stimvals.rf(2);
            Expt.Trials(id(j)).Op = 0;
            Expt.Trials(id(j)).Pp = 0;
            Expt.Trials(id(j)).st = 0;
        end
        end
    end
    
    if strcmp(Expt.Stimvals.et,'fx') && strcmp(Expt.Stimvals.e2,'fy')
        emname = strrep(Expt.Header.Name,'.mat','.eyecal.mat');
        if ~exist(emname,'file')
            MakeEyeCalfile(Expt,emname);
        else
            load(emname);
            fprintf('Eye gains RH %.2f LH %.2f RV %.2f LV %.2f\n',gains(1),gains(2),gains(3),gains(4));
        end
    end
%     if strmatch(et,'Dc')
%         ors = cat(2,Expt.Trials.or);
%     end
        Expts{nexpts} = Expt;
    nexpts = nexpts+1;
    end
%  max(spkids);
end
    
function  gains = MakeEyeCalfile(cExpt,emname);
 
emfile = strrep(cExpt.Header.Name,'.mat','.em.mat');
if ~exist(emfile,'file')
    gains = NaN;
    return;
end
load(emfile);
if ~isfield(Expt.Trials,'fx')
    for j = 1:length(cExpt.Trials)
        id = find([Expt.Trials.id] == cExpt.Trials(j).id);
        Expt.Trials(id).fx = cExpt.Trials(j).fx;
        Expt.Trials(id).fy = cExpt.Trials(j).fy;
    end
end
[gains, positions] = EyeCal(Expt,'ids',[cExpt.Trials.id]);
save(emname,'gains','positions');


function d = CreationDate(Text)

did = strmatch('uf',Text.text);
if isempty(did) %online file
    did = strmatch('bt',Text.text);
end
d = 0;
if length(did)
    for j = 1:length(did)
        ds = Text.text(did(j),:);
        dsid = strfind(ds,'Creat');
        if length(dsid)
            d = datenum(ds(dsid(1)+8:end));
            break;
        end
    end
end
if d == 0 %still not found
    did = strmatch('vebinoc',Text.text);
    if ~isempty(did)
        ds = Text.text(did(1),:);
        did = strfind(ds,' ');
         if length(did)
            d = datenum(ds(did(1)+1:end));
        end
    end
end

function [rfstr, rf] = MkUfl(name, Text, varargin)
%first make .ufl file with rf boxes, so that can build pen maps the
%old way
overwrite = 0;
j = 1;
rfstr = [];
rf = [];
while j <= length(varargin)
    if strncmpi(varargin{j},'overwrite',5)
        overwrite = 1;
    end
    j = j+1;
end
ufl = strrep(name,'.mat','.ufl');
rid = strmatch('cm=rf',Text.text);

AddTxtFile = strrep(name,'.mat','Add.txt');
fid = fopen(AddTxtFile,'r');
if fid > 0
    a = textscan(fid,'%d %s','delimiter','\n');
    fclose(fid);
    id = find(a{1} == -1);
    rfstrs = a{2}(id);
else
    rfstrs = {};
end

if isempty(rid) && isempty(rfstrs)
    return;
end
%a = textscan(Text.text(id,:),'cm=rf%f,%f:%fx%f,%fdeg pe%d %f %f%*s');
% trailing spaces seem to mess this up. text(id,1:65) works for most line
% but still barfs if a line is the wrong length
if isempty(rfstrs)
for j = 1:length(rid)
    a = sscanf(Text.text(rid(j),:)','cm=rf%f,%f:%fx%f,%fdeg pe%f %f,%f fx=%f,fy=%f');
    rfs(j,1:length(a)) = a;
end
else
    for j = 1:length(rfstrs)
        a = sscanf(rfstrs{j},'cm=rf%f,%f:%fx%f,%fdeg pe%f %f,%f fx=%f,fy=%f');
        rfs(j,1:length(a)) = a;
    end
end
% find lines suggesting RF was changed after a quantitative measure
oid = strmatch('RO',Text.text);
pid = strmatch('RP',Text.text);
sid = [oid pid strmatch('RO',Text.text)];
if length(sid)
    id = find(rid > max(sid))
end
for j = 1:size(rfs,2)
    rf(j) = mode(rfs(:,j));
end
if size(rfs,2) < 10
  rfstr = 'Missing RF data';
else
    rfstr = sprintf('cm=rf%.2f,%.2f:%.2fx%.2f,%.0fdeg pe%.0f %.1f,%.1f fx=%.2f,fy=%.2f\n',...
        rf(1),rf(2),rf(3),rf(4),rf(5),mode(rfs(:,6)),...
        mode(rfs(:,7)),mode(rfs(:,8)),mode(rfs(:,9)),mode(rfs(:,10)),rf);
end




if exist(ufl,'file') & ~overwrite
    return;
end


d = CreationDate(Text);
ds = [];
if d > 0
    ds = ['Created: ' datestr(d,'mm/dd/yyyy')];
end
of = fopen(ufl,'w');
if of > 0 
    fprintf(of,'%s\n',rfstr);
    for j = 1:length(rid)
%        fprintf(of,'%s\n',Text.text(id(j),:));
    end
    if ~isempty(ds)
        fprintf(of,'%s\n',ds);
    end
    fclose(of);
else
    questdlg(sprintf('Can''t Write %s',ufl),'test','OK','OK');
end

function Idx = AddError(err,Idx, show)

if ~isfield(Idx,'errs') | isempty(strmatch(err,Idx.errs)) 
    if show
        msgbox(err,'APlaySpkFile Error!!','modal');
    else
        fprintf([err '\n']);
    end
    Idx.errs = {Idx.errs{:} err};
    Idx.newerrs = Idx.newerrs+1;
else
    fprintf([err '\n']);
end

function FindMissingTimes(Events, Text, bstimes, estimes)

bid = strmatch('bss',Text.text);
bsstimes = Text.times(bid);
x = bsstimes-bstimes(1:length(bsstimes));
plot(bsstimes,x,'o');



function name = BuildName(name)
if isempty([strfind(name,'/') strfind(name,'\')])
    name = [pwd '/' name];
end
id = strfind(name,':');
if id
    name = name(id(1)+1:end);
else
    name = name;
end
if isunix
    name = strrep(name,'\','/');
end


