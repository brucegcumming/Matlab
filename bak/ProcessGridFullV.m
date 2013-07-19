
function DATA = ProcessGridFullV(name, varargin)
%DATA = ProcessGridFullV(name, varargin)
%make FullV files from .ns5 files. Apply clusters
%ProcessGridFullV(name,'online')
%ProcessGridFullV(name,'initial') =  'chopfile','BuildV','mains', 'submean', 'refcut', 'prebuild';
%
%ProcessGridFullV(name,...'parallel') uses parfor loops where possible
%(currently for refcuts only)

DATA = [];
Trials = [];
args = {'noerrs'};
gargs = {};
reindex = 0;
plotexpts = [];
expts = [];
probes = [];
BuildV = 0;
CheckV = 0;
BuildLFP = 0;
smoothw = 0;
submains = 0;
chopfile = 1;
gaplen = 40000;
autocut = 0;
fullVname = [];
buildmean = 0;
submean = 0;
savefile = 1;
sumv = [];
preperiod = 0;
Clusters = {};
refcut = 0;
prebuild = 0;
forcebuild = 0;
waitforfiles = 0;
highpass = NaN;
parallel = 0;
savespikes = 1;
allvargs = {};
buildargs = {};
j = 1;
while j <= length(varargin)
    if isfield(varargin{j},'Trials')
        Trials = varargin{j};
    elseif isfield(varargin{j},'ddelay')  %INdex file for nsx->expt
        DATA.idx = varargin{j};
    elseif isfield(varargin{j},'Events')  %INdex file for nsx->expt
        All = varargin{j};
    elseif iscell(varargin{j})
        Expts = varargin{j};
    elseif strncmpi(varargin{j},'autocut',5)
        autocut = 1;
    elseif strncmpi(varargin{j},'buildlfp',7)
        BuildLFP = 1;
    elseif strncmpi(varargin{j},'buildv',6)
        BuildV = 1;
    elseif strncmpi(varargin{j},'buildmean',6)
        buildmean = 1;
    elseif strncmpi(varargin{j},'buildnsmean',9)
        buildmean = 2;
    elseif strncmpi(varargin{j},'calcmains',9)
        submains = 2;
    elseif strncmpi(varargin{j},'checkpwr',6)
        CheckV = 1;
    elseif strncmpi(varargin{j},'expts',5)
        j = j+1;
        expts = varargin{j};
    elseif strncmpi(varargin{j},'chopfile',5)
        chopfile = 1;
    elseif strncmpi(varargin{j},'highpass',5)
        j = j+1;
        highpass = varargin{j};
    elseif strncmpi(varargin{j},'initial',5)
        chopfile = 1;
        BuildV = 1;
        submains = 1;
        submean = 1;
        refcut = 1;
        if regexp(name,'^[A-H]:') %local drive
            prebuild = 0;
        else
            prebuild = 1;
        end
    elseif strncmpi(varargin{j},'online',5)
        chopfile = 1;
        BuildV = 1;
        submains = 1;
        submean = 1;
    elseif strncmpi(varargin{j},'nochopfile',5)
        chopfile = 0;
    elseif strncmpi(varargin{j},'nosave',5)
        savespikes = 0;
    elseif strncmpi(varargin{j},'forcebuildv',8)
        BuildV = 1;
        forcebuild = 1;
        buildargs = {buildargs{:} varargin{j}};
    elseif strncmpi(varargin{j},'meanv',5)
        j = j+1;
        sumv = varargin{j};
    elseif strncmpi(varargin{j},'mains',5)
        submains = 1;
    elseif strncmpi(varargin{j},'maxgap',5)
        j = j+1;
        gaplen = varargin{j};
    elseif strncmpi(varargin{j},'nosave',5)
        savefile = 0;
    elseif strncmpi(varargin{j},'prebuildforce',10)
        prebuild = 2;
    elseif strncmpi(varargin{j},'parallel',4)
        parallel= 1;
    elseif strncmpi(varargin{j},'prebuild',5)
        prebuild = 1;
    elseif strncmpi(varargin{j},'probes',5)
        j = j+1;
        probes = varargin{j};
    elseif strncmpi(varargin{j},'plotexpts',5)
        j = j+1;
        plotexpts = varargin{j};
        gargs = {gargs{:} varargin{j-1} varargin{j}};
    elseif strncmpi(varargin{j},'smooth',5)
        j = j+1;
        smoothw = varargin{j};
    elseif strncmpi(varargin{j},'submean',5)
        submean = 1;
    elseif strncmpi(varargin{j},'refcut',5)
        refcut = 1;
    elseif strncmpi(varargin{j},'reindex',5)
        reindex = 1;
        gargs = {gargs{:} varargin{j}};
    elseif strncmpi(varargin{j},'relist',5)
        args = {args{:} varargin{j}};
    elseif strncmpi(varargin{j},'verbose',5)
        allvargs = {allvargs{:} varargin{j}};
    elseif strncmpi(varargin{j},'wait',4)
        waitforfiles = 1;
    elseif regexp(varargin{j},'Expt.*FullV.mat')
        fullVname = varargin{j};
    end
    j = j+1;
end

if submean && buildmean == 0
    buildmean == 3;
end
if isdir(name)
    datadir = name;
else
    datadir = fileparts(name);
end


ready = 0;
ts = 0;
while waitforfiles && ready == 0
    if ~exist(name)
        ready = 0;
    else
        ready = 1;
        d = dir([datadir '/*.nev']);
        for j = 1:length(d)
            nsfile = strrep([datadir '/' d(j).name],'.nev','.ns5');
            if ~exist(nsfile)
                if ready == 1 %first find
                    fprintf('Missing %s\n',nsfile);
                end
                ready = 0;
            end
        end
    end
    if ready == 0 
        fprintf('Waiting for files to be ready at %s\n',datestr(now));
        pause(600); %%10 mins
    end
end


if isempty(strfind(path,'BlackRock'))
    path(path,'/bgc/bgc/matlab/BlackRock');
end




if buildmean == 1
    for e = expts
    for j = 1:length(probes)
        p = probes(j);
        load(sprintf('%s/Expt%d.p%dFullV.mat',datadir,e,p));
        FullV.V = double(FullV.V);
        if j == 1
            sumv = FullV.V./std(FullV.V);
        else
            sumv(1:length(FullV.V)) = sumv(1:length(FullV.V)) + FullV.V./std(FullV.V);
        end
    end
    end
    sumv = sumv./length(probes);
    if savefile
        save(sprintf('%s/Expt%dFullVmean.mat',datadir,e),'sumv')
    end
    if submean
        X.idx.datdir = datadir;
        SubtractMean(X, e, probes, sumv,'verbose');
    end
    DATA = sumv;
    return;
end

if buildmean == 2
    ts = now;
    for e = expts
        id = find(DATA.idx.expt == e);
        filename = [DATA.idx.datdir '/' DATA.idx.names{id(1)}];
        filename = strrep(filename,'nev','ns5');
    for j = 1:length(probes)
        p = probes(j);
        tic
        v = openNSx('read',filename,['e:' num2str(p)],'p:int16');
        toc
        if j == 1
            sumv = double(v.Data);
        else
            sumv = sumv + double(v.Data);
        end
        vscale = v.MetaTags.VoltsPerDigit(p);
    end
    end
    fprintf('Took %.2fsec\n',mytoc(ts));
    sumv = sumv./length(probes);
    save(sprintf('%s/Expt%dFullVmean.mat',datadir,e),'sumv')
    DATA = sumv;
    return;
end


if ~isempty(Trials)
elseif isdir(name)  %online
    Trials.bstimes = [];
    Trials.mainstimes = [];
    Trials.Trials.Result = [];
    Trials.Trials.Start = [];
    Trials.Trials.End = [];
    Trials.DigMark = [];
    All.Events.codes = [];
    All.Events.times = [];
    nx = 0;
    d = dir([name '/Expt*.mat']);
    for j = 1:length(d)
        if isempty(strfind(d(j).name,'idx')) && isempty(strfind(d(j).name,'FullV')) ...
                && ~isempty(regexp(d(j).name,'Expt[0-9]*.mat'))
            [T, E, A] = APlaySpkFile([name '/' d(j).name], 'nospikes', args{:});
            nx = nx+1;
            Expts{nx} = E{1};
            Trials.DataType = T.DataType;
            Trials.Trials.Result = cat(2,Trials.Trials.Result,T.Trials.Result);
            Trials.Trials.Start = cat(2,Trials.Trials.Start,T.Trials.Start);
            Trials.Trials.End = cat(2,Trials.Trials.End,T.Trials.End);
            Trials.bstimes = cat(1,Trials.bstimes,T.bstimes);
            Trials.DigMark = cat(1,Trials.DigMark,T.DigMark);
            if isfield(T,'mainstimes')
                Trials.mainstimes = cat(1,Trials.mainstimes,T.mainstimes);
            end
            All.Events.codes = cat(1,All.Events.codes,A.Events.codes);
            All.Events.times = cat(1,All.Events.times,A.Events.times);
            Expts{nx}.loadname = [name '/' d(j).name];
        end
    end
    
else
    [Trials, Expts, All] = APlaySpkFile(name, 'nospikes',  args{:});
    for j = 1:length(Expts)
        Expts{j}.loadname = name;
    end
end


if ~isempty(fullVname)
    
    outname = [datadir '/' fullVname];
    if submains
        mcycle = [];
        if isempty(probes)
            load(outname);
            [FullV, mcycle] = RemoveMains(FullV, Trials.mainstimes,mcycle);
            if savefile
                fprintf('Saving %s at %s\n',outname,datestr(now));
                save(outname,'FullV');
            end
        else
        for p = probes
            name = regexprep(fullVname,'.p[0-9]*FullV',sprintf('.p%dFullV',p));
            load(outname);
            [FullV, mcycle] = RemoveMains(FullV, Trials.mainstimes,mcycle);
            if savefile
                fprintf('Saving %s at %s\n',outname,datestr(now));
                save(outname,'FullV');
            end
            DATA.mtrig(p,:) = FullV.mainsavg;
            plot(FullV.mainsavg);
            hold on;
        end
        end
    end
    return;
end

if isempty(Trials)
    DATA.idx = BuildGridIndex(name, [], gargs{:});
    expts = unique(DATA.idx.expt);
    if isfield(DATA.idx,'nprobes')
        probes = 1:max(DATA.idx.nprobes);
    else
        np = sscanf(Trials.DataType,'GridData %d');
    end
else
    np = sscanf(Trials.DataType,'GridData %d');
    if isempty(np)
        np = 96; % default
    end

    for j = 1:np
        DATA.probes(j).probe = j;
    end

    if isempty(probes)
        probes = [DATA.probes.probe];
    end
    if isempty(expts)
        expts = 1:length(Expts);
    else
        expts = expts(expts <= length(Expts));
    end

if chopfile
    preperiod = 10000;
    postperiod = 2000;
  nblk = 0;
  id = find(Trials.Trials.Result ~= 0);
  starts(1) = Trials.Trials.Start(id(1)) - preperiod;
  for j = 2:length(id)
      dt = Trials.Trials.Start(id(j))-Trials.Trials.End(id(j-1));
      if dt > gaplen
          nblk = nblk+1;
          ends(nblk) = Trials.Trials.End(id(j-1))+postperiod;
          starts(nblk+1) = Trials.Trials.Start(id(j))-preperiod;
      end
  end
  DATA.totaldur = Trials.Trials.End(end)+postperiod - starts(1);
  ends(nblk+1) = Trials.Trials.End(end)+postperiod;
  DATA.blockdur = sum(ends-starts)./10000;
  fprintf('%d blocks: %.1f/%.1f sec\n',nblk,DATA.blockdur,DATA.totaldur./10000);
  DATA.blocks = [starts ends];
  starts = starts./10000;
  ends = ends./10000;
end

sonid = find(All.Events.codes(:,1) ==48); %'0' = storage on
soffid = find(All.Events.codes(:,1) ==49); %'1' = storage off
for nexp = 1:length(Expts)
    if isfield(Trials,'bstimes')
        id = find(Trials.bstimes > Expts{nexp}.Header.Start & Trials.bstimes < Expts{nexp}.Header.End);
        Expts{nexp}.bstimes = Trials.bstimes(id);
    end
    if length(sonid)
        id = find(All.Events.times(sonid) > Expts{nexp}.Header.trange(1) & ...
            All.Events.times(sonid) < Expts{nexp}.Header.trange(2));
        Expts{nexp}.gridstoreon = All.Events.times(sonid(id));
        id = find(All.Events.times(soffid) > Expts{nexp}.Header.trange(1) & ...
            All.Events.times(soffid) < Expts{nexp}.Header.trange(2));
        Expts{nexp}.gridstoreoff = All.Events.times(soffid(id));
    end
    if ~isfield(Expts{nexp},'DigMark')
    if isfield(Trials,'DigMark')
        if Expts{nexp}.DigMark.codes(end) ~= 2 && ...
                Expts{nexp}.DigMark.times(end) < Expts{nexp}.Header.trange(2)./10000
            fprintf('Expt %d Doesn''t end with Stop\n',nexp);
        end
        
    elseif length(Expts{nexp}.gridstoreon)
        Expts{nexp}.DigMark = [];
        Expts{nexp}.DigMark.times = Expts{nexp}.gridstoreon;
        Expts{nexp}.DigMark.codes(1:length(Expts{nexp}.gridstoreon)) = 1;
    end
    end
end
if isfield(Trials,'DigMark')
    gargs = {'DigMark' Trials.DigMark gargs{:}};
end
DATA.idx = BuildGridIndex(name, Expts, gargs{:});
DATA.Expts = Expts;
end

    gotexpts = expts(ismember(expts,DATA.idx.expt));
    missedexpts = expts(~ismember(expts,DATA.idx.expt));
    for j = missedexpts
        fprintf('Missing Nsx/Nev files for Expt %d:%s %.1f\n',j,...
            Expts{j}.Header.expname,Expts{j}.Header.Start/10000);
    end
if prebuild
    for e = gotexpts
        id = find(DATA.idx.expt == e);
        filename = [DATA.idx.datdir '/' DATA.idx.names{id(1)}];
        filename = strrep(filename,'nev','ns5');
        needv = [];
        for k = probes
            rawfile = strrep(filename,'.ns5',['rawp' num2str(k) '.mat']);
            if ~exist(rawfile,'file')
                needv(k) = 1;
            end
        end
        if sum(needv) || prebuild == 2
            ts = now;
            nsx = openNSx('read',filename,'p:int16');
           DATA.readtime(e,1) = mytoc(ts);
            ts = now;
           
        for k = 1:size(nsx.Data,1)
            Data = nsx.Data(k,:);
            MetaTags = nsx.MetaTags;
            probe = k;
            rawfile = strrep(filename,'.ns5',['rawp' num2str(k) '.mat']);
            save(rawfile,'Data','MetaTags','probe');
        end
        clear v;
        clear V;
        DATA.readtime(e,2) = mytoc(ts);
        end
    end
end

DATA.prebuild = prebuild;
DATA.smoothw = smoothw;
DATA.chopfile = chopfile;
DATA.submean = submean;
if chopfile
DATA.starts = starts;
DATA.ends = ends;
DATA.gaplen = gaplen;
DATA.postperiod = postperiod;
end
DATA.preperiod = preperiod;

    Arrays = {};
    online = 0;
if BuildLFP
    for ex = 1:length(gotexpts)
        e = gotexpts(ex);
        id = find(DATA.idx.expt == e);
        tfirst = (Expts{e}.Trials(1).Start(1) - preperiod)/10000;
        outfile = [DATA.idx.datdir '/Expt' num2str(e) 'LFP.mat'];
        LFP.nsnames = strrep(DATA.idx.names{id},'nev','ns2');
        LFP.expname = Expts{e}.loadname;
        if isempty(xid)
            FullV.blklen = [];
            FullV.blkstart = [];
        end
        for k = 1:length(id)
            filename = [DATA.idx.datdir '/' DATA.idx.names{id(k)}];
            filename = strrep(filename,'nev','ns2');
            v = openNSx('read',filename,['e:' num2str(p)]);
            LFP.LFP = v.Data;
            LFP.tstart(k) = DATA.idx.toff(id(k))./10000;
        end
        save(outfile,'LFP');
    end
end

if BuildV
    DATA.prebuild = prebuild;
    DATA.submains = submains;
    DATA.buildmean = buildmean;
    DATA.chopfile = chopfile;
    fprintf('Building FullV for expts %s\n',sprintf('%d ',gotexpts));
    outvarname = 'FullV';
    if length(gotexpts) == 1  && online  %off for now
        e = gotexpts;
        DATA.tfirst = (Expts{e}.Trials(1).Start(1) - preperiod)/10000;
        if submains
            X = LoadFullV(DATA, Expts, e, 1);
            ts = now;
            [X, mcycle] = RemoveMains(X, Trials.mainstimes,[]);
            fprintf('First Mains subtration took %.2f\n',mytoc(ts));
            clear X;
        end
        for (p = probes)
            FullV = LoadFullV(DATA, Expts, e, p);
            fprintf('Loaded %s, p%d, says %d\n',FullV.name,p,FullV.chspk);
            if submains
                FullV = RemoveMains(FullV, Trials.mainstimes, mcycle);
            end
        end
        if buildmean
            if isempty(sumv)
                ts = now;
                sumv = ProcessGridFullV(DATA.idx.datdir,'expts',e,Trials,Expts,DATA.idx,'buildmean','probes',probes);
            end
            for p = probes
            end
        end
    else
            
        
%    parfor (ex = 1:length(gotexpts))
    for ex = 1:length(gotexpts)
        e = gotexpts(ex);
 %       BuildGridFullV(DATA, Expts, e, probes, buildargs{:});
        vid  = [];
        sumtrig = [];
        onetrig = [];
        id = find(DATA.idx.expt == e);
        tfirst = (Expts{e}.Trials(1).Start(1) - preperiod)/10000;
%        DATA.tfirst = tfirst; can't do this in parfor
        
        mcycle = [];
        xid = [];
        sumv = [];
        outvarname = 'FullV';
        if ~forcebuild
            doprobes = [];
            for (p = probes)
                outfile = [DATA.idx.datdir '/Expt' num2str(e) '.p' num2str(p) 'FullV.mat'];
                if ~exist(outfile)
                    doprobes = [doprobes p];
                end
            end
        else
            doprobes = probes;
        end
        for (p = probes)
            xid = [];
            FullV.V = [];
            ts = now;
            FullV.nsnames = DATA.idx.names(id);
            FullV.expname = Expts{e}.loadname;
            if isempty(xid)
                FullV.blklen = [];
                FullV.blkstart = [];
            end
            for k = 1:length(id)
                filename = [DATA.idx.datdir '/' DATA.idx.names{id(k)}];
                filename = strrep(filename,'nev','ns5');
                if prebuild
                   rawfile = strrep(filename,'.ns5',['rawp' num2str(p) '.mat']);
                    v = load(rawfile);
                    v.Data = double(v.Data);
                else
                v = openNSx('read',filename,['e:' num2str(p)]);
                end
                if smoothw > 0
                    sm = smooth(v.Data,smoothw);
                    v.Data = v.Data - sm;
                end
                if size(v.MetaTags.VoltsPerDigit,1) < p
                    vscale = v.MetaTags.VoltsPerDigit(1);
                    FullV.elecxy(1) = NaN;
                    FullV.elecxy(2) = NaN;
                else
                    vscale = v.MetaTags.VoltsPerDigit(p);
                    Array = GetArrayConfig(v.MetaTags);
                    FullV.elecxy(1) = Array.X(p);
                    FullV.elecxy(2) = Array.Y(p);
                end
                Arrays{e} = Array;
                FullV.V((1+length(FullV.V)):(length(v.Data)+length(FullV.V))) = v.Data;
                if isempty(xid)
                    FullV.blklen(k)= length(v.Data);
                    FullV.blkstart(k) = DATA.idx.toff(id(k))./10000;
                end
            end
            if size(v.MetaTags.ElecLabel,1) < p
                FullV.NSlabel = v.MetaTags.ElecLabel(1,:);
            else
                FullV.NSlabel = v.MetaTags.ElecLabel(p,:);
            end
            FullV.samper = 1./30000.237;
            if ~isempty(xid)
               FullV.V = FullV.V(vid);
            elseif chopfile && length(FullV.blkstart) ==1 && isempty(xid)
                tstart = FullV.blkstart;
                tlast = FullV.blkstart+FullV.blklen.*FullV.samper;
                if tfirst > tstart
                    firstsample = round((tfirst-FullV.blkstart)./FullV.samper);
                    xid = 1:firstsample-1;
                    FullV.blkstart = tfirst;
                    npre = 1;
                else
                    npre = 1;
                    firstsample = 2;
                end
               bid = find(starts > tfirst & starts > tstart & starts < tlast);
               for k = 1:length(bid)
                   b = round((starts(bid(k))-tstart)./FullV.samper);
                   a = round((ends(bid(k)-1)-tstart)./FullV.samper);
                   if (a> 0)
                       xid = [xid a:b];
                       if k+npre > 1
                           FullV.blklen(k+npre-1) = a-firstsample;
                       end
                   else
                       FullV.blklen(k+npre-1) = b-firstsample;
                   end
                   FullV.blkstart(k+npre) = starts(bid(k));
                   firstsample = b+1;
               end
               FullV.blklen(k+npre) = 1+length(v.Data)-firstsample;
               vid = setdiff(1:length(v.Data),xid);
               FullV.V = FullV.V(vid);
               FullV.chopratio = length(xid)./length(v.Data);
               FullV.preperiod = preperiod;
               FullV.postperiod = postperiod;
               FullV.chopgap = gaplen;
               lengths(p) = length(vid);
            end
            vm = max(abs(FullV.V.*vscale));
            FullV.intscale = [vm v.MetaTags.Resolution];
            outfile = [DATA.idx.datdir '/Expt' num2str(e) '.p' num2str(p) 'FullV.mat'];
            FullV.samper = 1./30000.237;
            FullV.name = outfile;
            FullV.builddate = now;
            FullV.chspk = p;
            FullV.exptno = e;
            FullV.buildtime = mytoc(ts);
            FullV.start = FullV.blkstart(1);
            if isfield(Expts{e},'loadname')
                FullV.matfile = Expts{e}.loadname;
            else
                FullV.matfile = name;
            end
            if submains == 2
                if isempty(mcycle)
                    [mainsfullv, mcycle] = RemoveMains(FullV, Trials.mainstimes,[],'calconly');
                end
            elseif submains %if 2, just calculated for first channel
                if isempty(mcycle)
                    [FullV, mcycle] = RemoveMains(FullV, Trials.mainstimes,[]);
                else
                    FullV = RemoveMains(FullV, Trials.mainstimes, mcycle);
                end
            end
            if buildmean == 3
                if isempty(sumv)
                    sumv = FullV.V./std(FullV.V);
                else
                    sumv(1:length(FullV.V)) = sumv(1:length(FullV.V)) + FullV.V./std(FullV.V);
                end
            end
                
            xscale = v.MetaTags.Resolution./max(abs(FullV.V));
            FullV.highpass = highpass;
            if isnan(highpass)
                [ratio, details] = HighLowRatio(FullV);
                if ratio > 10
                    FullV.highpass = 100;
                    fprintf('Low/High ratio %.2f - will need highpass\n',ratio);
                else
                    FullV.highpass = 0;
                    fprintf('Low/High ratio %.2f - No need for highpass\n',ratio);
                end
                FullV.HLratio = ratio;
                FullV.coilnoiseratio = details.coilratio;
            end
            FullV.V = int16(FullV.V * xscale);
            FullV.savetime = now;
            fprintf('Writing P%d(%d) %s (%s) (took %.2f at %s',p,FullV.chspk,outfile,outvarname,mytoc(ts),datestr(now));
            ts = FullV.savetime;
            SaveFullV(outfile, FullV);
            fprintf('+%.2f)\n',mytoc(ts));
        end
        if buildmean == 3
            if isempty(sumv)
                ts = now;
                sumv = ProcessGridFullV(DATA.idx.datdir,'expts',e,Trials,Expts,DATA.idx,'buildmean','probes',probes);
            end
            SubtractMean(DATA,e, probes, sumv);
            fprintf('Mean subtratction took %.1f at %s\n',mytoc(ts),datestr(now));
            comparenoise = 0;
            if comparenoise
            p = 1;
            outfile = [DATA.idx.datdir '/Expt' num2str(e) '.p' num2str(p) 'FullV.mat'];
            load(outfile);
            for j = length(mcycle):-1:1
                id = find(mcycle == j);
                sumtrig(j) = mean(sumv(id));
                onetrig(j) = mean(FullV.V(id));
            end
            plot(sumtrig,onetrig)
            end
        end
    end
end
end


id = [];
for j = 1:length(Arrays)
 if ~isempty(Arrays{j})
    id = j;
    break;
 end
end
if ~isempty(id)
    ArrayConfig = Arrays{id(1)};
    fprintf('Saving Array Config\n');
    save([DATA.idx.datdir '/ArrayConfig.mat'],'ArrayConfig');
else
    fprintf('Missing Array Config data\n');
end


if refcut
    if isempty(Clusters)
        cfile = [DATA.idx.datdir '/RefClusters.mat'];
        C = load(cfile);
        Clusters = C.Clusters;
        clear C;
        marked = CellToMat(Clusters,'marked');
        fprintf('%d Good cells, %d bad probes\n',sum(marked ==2),sum(marked ==3))
    end
    DATA.allvargs = {'matchcounts' 'nowatch' 'nocheck' 'noninteractive' allvargs{:}};
    if savespikes
        DATA.allvargs = {DATA.allvargs{:} 'savespikes'};
    end
    if parallel
        parfor (j = 1:length(gotexpts))
            fprintf('####Expt%d is Lab %d of %d\n',gotexpts(j),labindex,numlabs);
            MakeCut(DATA, Clusters, gotexpts(j),probes);
        end
    else
        for (j = 1:length(gotexpts))
            MakeCut(DATA, Clusters, gotexpts(j),probes);
        end
    end
end
if autocut
    parfor (j = 1:length(gotexpts))
        e = gotexpts(j);
        outfile = [DATA.idx.datdir '/Expt' num2str(e) '.p1FullV.mat'];
        AllVPcs(outfile,'tchan',probes,'GridData','savespikes','autocutall');
    end
end

if CheckV
    for j = 1:length(gotexpts)
        for k = 1:length(probes)
            outfile = [DATA.idx.datdir '/Expt' num2str(gotexpts(j)) '.p' num2str(probes(k)) 'FullV.mat'];
            load(outfile);
            [ratios(j,k), b] = HighLowRatio(FullV);
            coilratios(j,k) = b.coilratio;
            fprintf('%s: %.3f %.3f\n',outfile,ratios(j,k),b.coilratio);
        end
    end
    hist(ratios(:));
    DATA.ratios = ratios;
    DATA.coilratios = coilratios;
end

function SubtractMean(DATA, e, probes, sumv, varargin)

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'verbose',5)
        verbose = 1;
    end
    j = j+1;
end
    
    sumsq = sumv*sumv';
    for p = probes
        outfile = [DATA.idx.datdir '/Expt' num2str(e) '.p' num2str(p) 'FullV.mat'];
        load(outfile);
        if ~isfield(FullV,'chspk')
            FullV.chspk= p;
        end
        if verbose
        fprintf('Writing P%d(%d) %s  at %s\n',p,FullV.chspk,outfile,datestr(now));
        end
        FullV.V = double(FullV.V);
        FullV.sumscale = (sumv(1:length(FullV.V))*FullV.V')./sumsq;
        FullV.V = FullV.V - sumv(1:length(FullV.V)) .* FullV.sumscale;
        FullV.V = int16(FullV.V);
        FullV.savetime(2) = now;
        SaveFullV(outfile, FullV);
    end



function MakeCut(DATA, Clusters, e, probes)

tag = sprintf('Expt%d',e);
for p = probes;
    outfile = [DATA.idx.datdir '/Expt' sprintf('%d.p%dFullV.mat',e,p)];
    if isfield(Clusters{p},'shape')
        if isfield(Clusters{p},'marked') && Clusters{p}.marked == 2
            res = AllVPcs(outfile,'tchan',p,'GridData','reapply',Clusters{p},DATA.allvargs{:},'toptag',tag);
        else
            res = AllVPcs(outfile,'tchan',p,'GridData','reapply',Clusters{p},DATA.allvargs{:},'toptag',tag);
        end
    else
        fprintf('Ref Cluster P%d is empty\n',p);
    end
end
if res.toplevel
close(res.toplevel);
end

function took = SaveFullV(name, FullV)

tic;
save(name,'FullV');
took = toc;


function [ratio, details] = HighLowRatio(FullV)
w = FullV.blklen(1);
if w > length(FullV.V)
    fprintf('blklen > length in %s\n,',FullV.name);
    w = length(FullV.V);
end
    x = abs(fft(double(FullV.V(1:w))));
    freq = 30000./w;
    freqs = freq * [1:w];
    lid = find(freqs > 1 & freqs < 10);
    hid = find(freqs > 1000 & freqs < 7000);
    cid = find(freqs > 5000 & freqs < 7000);
    [a,b] = max(x(cid)); %coil noise
    if b < 500
        hid = setdiff(hid,cid(1:b+500));
    elseif b > length(cid)-500
        hid = setdiff(hid,cid(b-500:end));
    else
        hid = setdiff(hid,cid(b-500:b+500));
    end
    
    ratio = mean(x(lid))/mean(x(hid));
    hid = find(freqs > 1000 & freqs < 5500);
    details.pwr = x;
    details.freq = freqs;
    details.coilratio = a./mean(x(hid));
    
    
function [FullV, MainsTrig] = RemoveMains(FullV, mt, MainsTrig, varargin)


calconly = 0;
newcycle = 0;
mcycle = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'calconly',8)
        calconly = 1;
    end
    j = j+1;
end

    
mp = prctile(diff(mt), 90)./10000;  %period of mains cycle in sec
ns = ceil(mp./FullV.samper); %# of samples per mains cycle.
mt = mt./10000;  %work in sec
if ~isfield(FullV,'t')
    first = 1;
    vt(sum(length(FullV.blklen))) = 0;
    for j = 1:length(FullV.blklen)
        last = first+FullV.blklen(j)-1;
        vt(first:last) = FullV.blkstart(j)+[1:FullV.blklen(j)].*FullV.samper;
        first = last+1;
    end
end

tst = now;
if isempty(MainsTrig)
    starts = [1 cumsum(FullV.blklen)];
    ends = starts(2:end)-1;
    mcycle(1:length(FullV.V)) = NaN;
    for j = 1:length(FullV.blklen)
        ts = vt(starts(j));
        te = vt(ends(j));
        ti = vt(starts(j):ends(j));
        mid = find(mt > ts & mt < te+0.02);
        k = 1;
        tic;
        nextt = mt(mid(k));
        if mid(k) == 1
            lastt = mt(mid(k))-0.0167;
        else
            lastt = mt(mid(k)-1);
        end
        for t = 1:length(ti)
            if ti(t) > nextt & k < length(mid)
                k = k+1;
                lastt = nextt;
                nextt = mt(mid(k));
            end
            mcycle(starts(j)+t-1) = ti(t)-lastt;
        end
    end
    FullV.mainsbuild(1) = mytoc(tst);
    mcycle = ceil(mcycle.*30000);
    mcycle = mcycle(1:length(FullV.V));
    MainsTrig.mcycle = mcycle;
    newcycle = 1;
end

if calconly
    for j = ns:-1:1
        id = find(mcycle == j);
        mtrig(j) = mean(FullV.V(id));
        mainsfullv(id) = mtrig(j);
    end
    FullV.fullv = mainsfullv;
    FullV.mtrig = mtrig;
    return;
end

adjgain = 1;
if adjgain
mmean = zeros(size(FullV.V));
rawv = FullV.V;
end
for j = ns:-1:1
    if newcycle
        MainsTrig.mids{j} = find(mcycle == j);
    end
    id = MainsTrig.mids{j};
    id = id(id < length(FullV.V));
    mtrig(j) = mean(FullV.V(id));
    mmean(id) = mtrig(j);
    FullV.V(id) = FullV.V(id)-mtrig(j);
end
FullV.mainsbuild(2) = mytoc(tst);
FullV.mainsavg = mtrig; 

function [FullV, mcycle] = RemoveMainsB(FullV, mt, mcycle, varargin)


calconly = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'calconly',8)
        calconly = 1;
    end
    j = j+1;
end
    
mp = prctile(diff(mt), 90)./10000;  %period of mains cycle in sec
mper = mp./FullV.samper; %# of samples per mains cycle - including fraction.
ns = ceil(mper);
mt = mt./10000;  %work in sec
if ~isfield(FullV,'t')
    first = 1;
    vt(sum(length(FullV.blklen))) = 0;
    for j = 1:length(FullV.blklen)
        last = first+FullV.blklen(j)-1;
        pts = first:last;

        vt(pts) = FullV.blkstart(j)+[1:FullV.blklen(j)].*FullV.samper;
        mid = find(mt < vt(first));
        a = mper*((vt(first)-mt(mid(end)))/mp);
        a = round(mod(a,mper));
        mcycle(pts) = round(mod(a+pts,mper));
        starts(j) = first;
        ends(j) = last;
        first = last+1;
    end
end

ts = now;

if calconly
    for j = ns:-1:1
        id = find(mcycle == j);
        mtrig(j) = mean(FullV.V(id));
        mainsfullv(id) = mtrig(j);
    end
    FullV.fullv = mainsfullv;
    FullV.mtrig = mtrig;
    return;
end

for k = 0:100
    offset = k*300000;
    M = mcycle(offset+[1:300000]);
    for j = ns:-1:1
        id = find(M == j)+offset;
        mtrig(j) = mean(FullV.V(id));
        FullV.V(id) = FullV.V(id)-mtrig(j);
    end
    mtrigs(k+1,:) = mtrig;

end
FullV.mainsbuild(2) = mytoc(ts);
FullV.mainsavg = mtrig; 



function OldMcycle()

for k = 1:length(mid)
        [a,b] = min(abs(ti-mt(mid(k))));
        ii = b+starts(j)-1;
        mcycle(ii:ii+ns+3) = 1:ns+4;
        if k == 1
            pre = ii-ns:ii-1;
            pre = pre(pre > 0);
            id = find(isnan(mcycle(pre)));
            mcycle(id) = ns-length(id)+1:ns;
        end
end

    

function Array = GetArrayConfig(MetaTags)
for j = 1:size(MetaTags.ElecLabel,1); 
    E(j) = sscanf(MetaTags.ElecLabel(j,:),'elec%d'); 
end
Array.Y = 1+mod(E-1,10);
Array.X = ceil(E/10);


function FullV = LoadFullV(DATA, Expts, e, p)

FullV.blklen = [];
FullV.blkstart = [];
FullV.V = [];
ts = now;
FullV.expname = Expts{e}.loadname;
tfirst = DATA.tfirst;
starts = DATA.starts;
id = find(DATA.idx.expt == e);
FullV.nsnames = DATA.idx.names(id);
    for k = 1:length(id)
                filename = [DATA.idx.datdir '/' DATA.idx.names{id(k)}];
                filename = strrep(filename,'nev','ns5');
                if DATA.prebuild
                   rawfile = strrep(filename,'.ns5',['rawp' num2str(p) '.mat']);
                    v = load(rawfile);
                    v.Data = double(v.Data);
                else
                v = openNSx('read',filename,['e:' num2str(p)]);
                end
                if DATA.smoothw > 0
                    sm = smooth(v.Data,DATA.smoothw);
                    v.Data = v.Data - sm;
                end
                vscale = v.MetaTags.VoltsPerDigit(p);
                Array = GetArrayConfig(v.MetaTags);
                FullV.elecxy(1) = Array.X(p);
                FullV.elecxy(2) = Array.Y(p);
                Arrays{e} = Array;
                FullV.V((1+length(FullV.V)):(length(v.Data)+length(FullV.V))) = v.Data;
                FullV.blklen(k)= length(v.Data);
                FullV.blkstart(k) = DATA.idx.toff(id(k))./10000;
            end
            FullV.NSlabel = v.MetaTags.ElecLabel(p,:);
            FullV.samper = 1./30000.237;
            if DATA.chopfile && length(FullV.blkstart) ==1
                tstart = FullV.blkstart;
                tlast = FullV.blkstart+FullV.blklen.*FullV.samper;
                if tfirst > tstart
                    firstsample = round((tfirst-FullV.blkstart)./FullV.samper);
                    xid = 1:firstsample-1;
                    FullV.blkstart = tfirst;
                    npre = 1;
                else
                    npre = 1;
                    firstsample = 2;
                end
               bid = find(starts > tfirst & starts > tstart & starts < tlast);
               for k = 1:length(bid)
                   b = round((starts(bid(k))-tstart)./FullV.samper);
                   a = round((DATA.ends(bid(k)-1)-tstart)./FullV.samper);
                   xid = [xid a:b];
                   FullV.blkstart(k+npre) = starts(bid(k));
                   if k+npre > 1
                   FullV.blklen(k+npre-1) = a-firstsample;
                   end
                   firstsample = b+1;
               end
               FullV.blklen(k+npre) = 1+length(v.Data)-firstsample;
               vid = setdiff(1:length(v.Data),xid);
               FullV.V = FullV.V(vid);
               FullV.chopratio = length(xid)./length(v.Data);
               FullV.preperiod = DATA.preperiod;
               FullV.postperiod = DATA.postperiod;
               FullV.chopgap = DATA.gaplen;
            end
            vm = max(abs(FullV.V.*vscale));
            FullV.intscale = [vm v.MetaTags.Resolution];
            outfile = [DATA.idx.datdir '/Expt' num2str(e) '.p' num2str(p) 'FullV.mat'];
            FullV.samper = 1./30000.237;
            FullV.name = outfile;
            FullV.builddate = now;
            FullV.chspk = p;
            FullV.exptno = e;
            FullV.buildtime = mytoc(ts);
            FullV.start = FullV.blkstart(1);
            if isfield(Expts{e},'loadname')
                FullV.matfile = Expts{e}.loadname;
            else
                FullV.matfile = name;
            end
