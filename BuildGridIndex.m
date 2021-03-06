function idx =BuildGridIndex(name, Expts, varargin)
%BUILDGRIDINDEX idx = BuildGridIndex(name, Expts, ...)
%idx = BuildGridIndex(name, Expts, ...)
%name is a filename or direcotry where the .mat file lives (used to
%construct path for where .nev files live
% Expts is a cell array of Expts, as returned by APlaySpkFile
% BuildGridIndex(filename, [], 'reindex') to rescan nev files
%
%Version 2 uses intervals between stimon/off pulse to match ns5 files with
%Expts.
% BuildGridIndex(filename, [], 'compare') Compares old vs new list 
% BuildGridIndex(filename, [], 'update') checks to see if the existing
% index was built the old way. If so, it makes a copy of the old file, rebuilds and compares the two.
% The returned structure has a filed oldidx with the oldx indices
% BuildGridIndex(filename, [], 'forcematches','jbeG086005.ns5',5, 1)  Forces the mathing expt for an ns5 file, with a 1 event offset


    idx = [];
   datdir = 'F:/Utah/jbe/';
   reindex = 0;
   nevdir = [];
   plotexpts = [];
   checktag = 'TrialCheck';
   preperiod = 5000;
   postperiod = 5000;
   compareidxs = 0;
   maxgap = 10000; %default is to chop if gap > 1 sec
   DigMark = [];
   plottype = 0;
   Trials = [];
   VERSION = 2.0;
   forcematches = {};
   forceexpts = [];
   forcebsoff = [];
   showerr = 1;
   aargs = {};
   usealltrials = 0;
if nargin == 1
    Expts = [];
end
   
j = 1;
while j <= length(varargin)
    if isfield(varargin{j},'ExptList')
        Trials = varargin{j};
    elseif strncmpi(varargin{j},'compare',5)
        compareidxs = 1;
    elseif strncmpi(varargin{j},'digmark',5)
        j = j+1;
        DigMark = varargin{j};
    elseif strncmpi(varargin{j},'forcematches',8)
        j = j+1;
        forcematches = {forcematches{:} varargin{j}};
        j = j+1;
        forceexpts = [forceexpts varargin{j}];
        j = j+1;
        forcebsoff = [forcebsoff varargin{j}];
    elseif strncmpi(varargin{j},'nevdir',5)
        j = j+1;
        nevdir = varargin{j};
    elseif strncmpi(varargin{j},'plotfiles',5)
        plottype = 1;
    elseif strncmpi(varargin{j},'plotexpts',5)
        j = j+1;
        plotexpts = varargin{j};
    elseif strncmpi(varargin{j},'noerrs',5)
        showerr = 0;
    elseif strncmpi(varargin{j},'plotidx',5)
        plottype = 3;
    elseif strncmpi(varargin{j},'reindex',5)
        reindex = 1;
    elseif strncmpi(varargin{j},'usealltr',8) %remake the new way (V.2)
        aargs= {aargs{:} varargin{j}};
        usealltrials =1;
    elseif strncmpi(varargin{j},'update',5) %remake the new way (V.2)
        reindex = 2;        
    end
    j = j+1;
end


if isstruct(name)
    PlotIdx(name);
    return;
elseif isdir(name)
       datdir = name;
   else
       datdir = fileparts(name);
end
if isempty(nevdir)
    nevdir = datdir;
end

   plotsummary = 2;

   if isempty(strfind(path,'BlackRock'))
      addpath([GetFilePath('bgcmatlab') '/BlackRock']);
   end
   idxfile = [datdir '/FileIdx.mat'];
   
   if ischar(Expts)
   end
   
   if exist(idxfile) && (compareidxs || reindex ==2)
       oldidx = load(idxfile);
       if isfield(oldidx.idx,'version') && oldidx.idx.version > 1
           fprintf('%s is already Version %.1f\n',oldidx.idx.version);
           idx = oldidx.idx;
           return;
       elseif reindex == 2 && ~isfield(oldidx.idx,'version')
           oldname = strrep(idxfile,'FileIdx','OldIdx');
           idx = oldidx.idx;
           save(oldname,'idx');
           reindex = 1;
           clear idx;
       end
       idx.oldidx = oldidx.idx;
   elseif exist(idxfile,'file') & ~reindex
       a = load(idxfile);
       idx = a.idx;
       idx.datdir = datdir;
       idx.nevdir = nevdir;
       nidx = length(idx.names);
       if plottype == 1
           PlotExptFiles(idx, Expts, Trials);
       end
       return;
   else
       nidx = 0;
   end
   
   if exist(idxfile,'file') && reindex
       BackupFile(idxfile,'print');
   end
   if isempty(Expts)
       if isdir(name) %online files
           [Expts, Trials] = ReadExptDir(name, 'online', 'noerrs');
           if iscell(Trials)
               for nex = 1:length(Expts)
                   Expts{nex}.bstimes = Trials{nex}.Trials.bstimes';
                   Expts{nex}.gridstoreon = Trials{nex}.Trials.bstimes(1);
               end
           end
           idx = BuildMatFiles(nevdir);
       else
           Expts = GetExpts(name, aargs{:});
       end
   end
   
   
   offtimes = [];
   ontimes = [];
   sampleoffs = [];
   sampleons = [];
   offids = [];
   cuts = [];
   ncut = 0;
   
   STOREBIT=2; %use 3 for older versions with no pausing

   for j = 1:length(Expts)
       starts(j) = Expts{j}.Header.CreationDate + Expts{j}.Header.Start./(10000 * 60 *60 *24);
       ebstimes{j} = Expts{j}.bstimes;
       id = find(Expts{j}.DigMark.codes ==1);
       if ~isfield(Expts{j},'gridstoreon') || isempty(Expts{j}.gridstoreon)
           Expts{j}.gridstoreon = Expts{j}.bstimes(1);;
       end
       if isempty(id)
       ontimes = [ontimes Expts{j}.gridstoreon./10000];
       else
           %are digmark times in sec or 1/10 msec? Online may differ from
           %final
       ontimes = [ontimes Expts{j}.DigMark.times(id(1))];
       end
       if id == 1
           %Online 2 files -> 2 Digmark strucutres. 
           if size(DigMark) >= j
               dj = j;
           else
               dj = 1;
           end
           sampleon = Expts{j}.Header.CreationDate + Expts{j}.DigMark.times(id)'./(60 *60 *24);
%           sampleon = Expts{j}.DigMark.times(id);
           if ~isempty(DigMark)
               tid = find(DigMark(dj).times == Expts{j}.DigMark.times(id));
               if tid > 1
                   sampleon = Expts{j}.Header.CreationDate + DigMark(dj).times(tid-1)./(60*60*24);
               end
           end
       elseif isempty(id)
           sampleon = Expts{j}.Header.CreationDate + Expts{j}.gridstoreon./(60 *60 *24*10000);
       else
           sampleon = Expts{j}.Header.CreationDate + Expts{j}.DigMark.times(id)'./(60 *60 *24);
       end
       if isempty(sampleon)
           sampleon = Expts{j}.Header.CreationDate + Expts{j}.bstimes(1)./(60 *60 *24*10000);
       end
       nson = length(sampleon);
       sampleons = [sampleons sampleon(1)];
       exptids(j) = Expts{j}.Header.idrange(end);
       id = find(Expts{j}.DigMark.codes ==2);
       
 %Nov 2012 Onwards, use Trial Isis to do matching, so start/end times are less critical
 %
       if isempty(id)
           fprintf('Expt %d Missing Sample off marker\n',j);
           if isempty(Expts{j}.gridstoreoff)
               Expts{j}.gridstoreoff(1) = Expts{j}.Header.End;
           end
           sampleoff = Expts{j}.Header.CreationDate + Expts{j}.gridstoreoff(1)./(60 *60 *24 *10000);
           offtimes = [offtimes Expts{j}.gridstoreoff(1)];
       else
           sampleoff = Expts{j}.Header.CreationDate + Expts{j}.DigMark.times(id)'./(60 *60 *24);
           offtimes = [offtimes Expts{j}.DigMark.times(id(1))'];
       end
       if length(id) > length(nson)
           sampleoffs = [sampleoffs sampleoff(1:length(nson))];
       else
           sampleoffs = [sampleoffs sampleoff];
       end
       offids = [offids j];
%       eestimes{j} = Expts{j}.estimes;
       for t = 2:length(Expts{j}.Trials)
           gaps(t) = Expts{j}.Trials(t).Start(1)-Expts{j}.Trials(t-1).End(end);
           if  gaps(t) > maxgap
            ncut = ncut+1;   
            cuts(ncut,1) = Expts{j}.Trials(t-1).End(end)+postperiod;
            cuts(ncut,2) = Expts{j}.Trials(t).Start(1)-preperiod;
           end
       end
   end
   ends(1:j-1) = starts(2:end);
   ends(j) = starts(j) + (max(ebstimes{j})./(10000 * 60 *60 *24));
   
   
%First Build a list of Nev files and their times
%This works for multiple files per expt. 
   d = dir([nevdir]);
   filenames = {d.name};
   %d = dir([datdir '/*.nev']);
   newf = 0;
   nnev = 0;
   usenev = 1;
   for j = 1:length(d)
       if ~isempty(strfind(d(j).name,'.nev')) && usenev %this has starttime and Dig Events
           nevfile = [nevdir '/' d(j).name];
           matfile = strrep(d(j).name,'.nev','.mat');
           mid = strmatch(matfile,filenames);
           if length(mid)
               agediff = d(j).datenum > d(mid).datenum;
           else
               agediff = 1;
           end
           if reindex || ~isfield(idx,'names') || isempty(strmatch(d(j).name,idx.names)) || ...
                   agediff > 0
               if  agediff > 0
                   nev = openNEV('read','nomat','noparse','nowarning',nevfile);
               else
                   nev = openNEV('read','noparse','nowarning','nowaves',nevfile);
               end
               eb = int16(bitand(3,nev.Data.SerialDigitalIO.UnparsedData));
               onoff = diff(int16(bitand(4,nev.Data.SerialDigitalIO.UnparsedData)));
               nnev =nnev+1;
               if isempty(nev.Data.SerialDigitalIO.TimeStampSec)
                   ns = 0;
                   bstimes = [];
                   estimes = [];
               else
                   id = find(onoff > 0);
                   bstimes = nev.Data.SerialDigitalIO.TimeStampSec(id+1).*10000;
                   id = find(onoff < 0);
                   estimes = nev.Data.SerialDigitalIO.TimeStampSec(id+1).*10000;
                   if diff(size(bstimes)) > 1
                       bstimes = bstimes';
                       estimes = estimes';
                   end
               end
                   ns = length(bstimes);
               ts = nev.MetaTags.DateTimeRaw;
               tstart  = datenum(ts(1),ts(2),ts(4),ts(5),ts(6),ts(7));
%When storage is turned off the lowest bit is dropped 1ms boefore bit2, 
% so lowest two bits == 2. Only other time this happens is at the stare
% when bit 1 is toggled to mark times.
               if isempty(nev.Data.SerialDigitalIO.TimeStampSec)
                   fprintf('%s Emptyt',nevfile);
                   requesttime = NaN;
                   tstop = NaN;
               elseif length(nev.Data.SerialDigitalIO.UnparsedData) == 0 || ...
                       bitand(3, nev.Data.SerialDigitalIO.UnparsedData(end)) > 0
                   fprintf('%s Missing Final Off event',nevfile);
                   [a,c] = min(abs(tstart - sampleons));
                   b = offids(c);
                   cpuclockdiff = (tstart-sampleons(c)).*(60*60*24);
                   nevoff =  nev.Data.SerialDigitalIO.TimeStampSec(end);
                   tstop = tstart+nevoff ./(24 * 60 * 60);
                   nevon =  nev.Data.SerialDigitalIO.TimeStampSec(1);
                   tdiff = ontimes(c)-nevon; %difference in sec
                   cpuclockdiff = (tstart-sampleons(c)).*(60*60*24);
                   requesttime = sampleons(c) * 10000;
               else
                   id = find(bitand(3,nev.Data.SerialDigitalIO.UnparsedData) ==2);
                   if isempty(id)
                       nevoff = 0;
                   else
                       nevoff =  nev.Data.SerialDigitalIO.TimeStampSec(id(end));
                   end
                   if nevoff < 1
                       fprintf('Nominally short file %.2f\n', nevoff);
                       nevoff =  nev.Data.SerialDigitalIO.TimeStampSec(end);
                       if isfield(nev.Data.Spikes,'TimeStamp') && ~isempty(nev.Data.Spikes.TimeStamp)
                           lastspk = double(nev.Data.Spikes.TimeStamp(end))./nev.MetaTags.TimeRes;
                           nevoff = max([nevoff lastspk]);
                       end
                       t = tstart+nevoff./(60 * 60 * 24);
                       xid = find(sampleoffs) > t;
                       if ~isempty(id)
                       fprintf('Missing ~ %.2f sec of data\n',sampleoffs(id(1))-1);
                       end
                   end
                   t = tstart+nevoff./(60 * 60 * 24);
                   tstop = t;
                   [a,c] = min(abs(t - sampleoffs));
                   b = offids(c);
                   tdiff = offtimes(c)-nevoff; %difference in sec
                   cpuclockdiff = (tstop-sampleoffs(c)).*(60*60*24);
                   requesttime = sampleons(c) * 10000;
               end
               
               dtid = 1:length(bstimes);
               dx = diff(bstimes);
               for e = 1:length(ebstimes)
                   if length(ebstimes{e}) >= length(bstimes) && length(bstimes) > 1
                       xsc = [];
                       for dt = 0:length(ebstimes{e})-length(bstimes)
                           xsc(dt+1) = std(diff(ebstimes{e}(dt+dtid))-dx);
                       end
                       [mindt(e,nnev), bsoffset(e,nnev)] = min(xsc);
                       if length(xsc) > 1
                           xsc = sort(xsc);
                           amindt(e,nnev) = xsc(2);
                       else
                           amindt(e,nnev) = NaN;
                       end
                   elseif length(ebstimes{e}) == length(bstimes)-1 && length(bstimes) > 1
                       xsc = std(diff(ebstimes{e})-dx(2:end)); %can get one extra in nev...
                       mindt(e,nnev) = xsc;
                       bsoffset(e,nnev) = -1;
                       amindt(e,nnev) = NaN;
                   elseif length(bstimes) > 1 %if more trials in nev that expt, can't be  a match.  
%But if some trials have been deleted from Expt (e.g. because of missing RC sequnece, can get this
%So add calculation.  But shouldn't need this, so for now, set mindt to Nan
                       xid = 1:length(ebstimes(e));
                       xsc = [];
                       for dt = 0:length(bstimes)-length(ebstimes{e});
                           xsc(dt+1) = std(diff(ebstimes{e}-dx(xid+dt)));
                       end
                       [mindt(e,nnev), bsoffset(e,nnev)] = min(xsc);
                       bsoffset(e, nnev) = bsoffset(e, nnev) * -1;
                       if length(xsc) > 1
                           xsc = sort(xsc);
                           amindt(e,nnev) = xsc(2);
                       else
                           amindt(e,nnev) = NaN;
                       end
                       mindt(e,nnev) = NaN;
                   else
                       bsoffset(e, nnev) = NaN;
                       amindt(e,nnev) = NaN;
                       mindt(e,nnev) = NaN;
                   end
               end
               nidx = nnev;
               [bestdt, b] = min(mindt(:,nnev));
               nextstd = min(amindt(:,nnev)); %lowest of second best matches
               if bestdt < 10 && (mindt(b,nnev)/nextstd < 0.001 || isnan(nextstd))
                   bsoff = bsoffset(b,nnev);
                   newf = newf+1;
                   idx.names{nidx} = d(j).name;
                   idx.expt(nidx) = b;
                   if bsoff < 0
                       idx.toff(nidx) = ebstimes{b}(1)-bstimes(1-bsoff);
                   else
                       idx.toff(nidx) = ebstimes{b}(bsoff)-bstimes(1);
                   end
                   idx.bsstd(nidx) = bestdt;
                   nt = length(bstimes);
                   idx.firstbs(nidx) = bsoff;
                   idx.start(nidx) = tstart + bstimes(1)./(10000 * 60 * 24); %datenum
                   idx.end(nidx) = tstart + estimes(end)./ (10000 * 60 * 24);
                   idx.nt(nidx) = length(estimes);
                   idx.bstimes{nidx} = bstimes;
                   idx.evtimes{nidx} = nev.Data.SerialDigitalIO.TimeStampSec .* 10000;
                   idx.digin{nidx} = bitand(7,nev.Data.SerialDigitalIO.UnparsedData);
                   cid = find(cuts(:,1) > idx.toff(nidx) & cuts(:,1) < idx.toff(nidx)+nevoff*10000);
                   idx.cuts{nidx} = cuts(cid,:);
                   needed(j) = 1;
               else
                   fprintf('Best Tscatter for %s %.2f, next best %,2f\n',d(j).name,mindt(b,nnev),nextstd)
                   idx.expt(nidx) = 0;
                   idx.tdiff(nidx) = tdiff;
                   idx.bsstd(nidx) = 0;
                   if ns
                   idx.digdt(nidx) = bsoffset(b,nnev);
                   else
                   idx.digdt(nidx) = NaN;
                   end
                   needed(j) = 0;
               end
               id = find(strcmp(d(j).name,forcematches))
               if length(id) ==1
                   b = forceexpts(id);
                   idx.expt(nidx) = b;
                   bsoff = forcebsoff(id);
                   figure;
                   if bsoff < 0
                       idx.toff(nidx) = ebstimes{b}(1)-bstimes(1-bsoff);
                       plot(diff(ebstimes{b}));
                       hold on;
                       plot(diff(bstimes(1-bsoff:end)),'r');
                   else
                       idx.toff(nidx) = ebstimes{b}(bsoff)-bstimes(1);
                       plot(diff(ebstimes{b}(bsoff:end)));
                       hold on;
                       plot(diff(bstimes),'r');
                   end
                   title(sprintf('Expt %d and %s, bso offset %d',b,d(j).name,bsoff));
               end
               idx.stds(nidx,1) = bestdt;
               idx.stds(nidx,2) = nextstd;
               if length(nev.Data.SerialDigitalIO.TimeStampSec)
                   idx.ddelay(nidx) = nev.Data.SerialDigitalIO.TimeStampSec(1).*10000;
               else
                   idx.ddelay(nidx) = NaN;
               end
               if length(nev.Data.SerialDigitalIO.TimeStampSec)
                   idx.ddelay(nidx) = nev.Data.SerialDigitalIO.TimeStampSec(1).*10000;
               else
                   idx.ddelay(nidx) = NaN;
               end
               if requesttime > 0 && c > 1
                   idx.offinterval(nidx) = sampleons(c)-offtimes(c-1);
               end
               idx.starttime(nidx) = tstart;
               idx.stoptime(nidx) = tstop;
               if ~isempty(nev.Data.SerialDigitalIO.TimeStampSec)
                   idx.lastevtime(nidx) =nev.Data.SerialDigitalIO.TimeStampSec(end);
               else
                   idx.lastevtime(nidx) =  NaN;
               end
               idx.names{nidx} = d(j).name;
               idx.cpuclockdiff(nidx) = cpuclockdiff;
               idx.bsoff(nidx) = bsoffset(b,nnev);
               if length(nev.Data.SerialDigitalIO.UnparsedData) 
                   idx.lastdio(nidx) = bitand(7,nev.Data.SerialDigitalIO.UnparsedData(end));
               end
           end
       end
   end
   idx.exbstimes = ebstimes;
   idx.version = VERSION;
   idx.builddate = now;
   idx.usealltrials = usealltrials;
   
   if  nidx == 0
       fprintf('Missing NeV Data Filesin %s\n',name);
       return;
   end
   
   missing = setdiff(1:length(Expts),idx.expt);
   if length(missing) 
       idx.missing = missing;
       idx = AddError(idx, sprintf('Missing NEV for Expts %s',sprintf('%d ',missing)));
       if showerr
           warndlg(idx.errs{end},'Expts Missing');
       end
   end
   eid = find(idx.expt > 0);
   badsd = find(idx.bsstd(eid) > 100);
   if length(badsd)
       idx.badsd = badsd;
       warndlg(sprintf(' Expts %s',sprintf('%d ',badsd)),'StimOn Markers Poor Match');
       for j = 1:length(badsd)
           iid = eid(badsd(j))
           fprintf('Expt %d offsets SD = %.1f\n',idx.expt(iid),idx.bsstd(iid));
       end
   end
   
idx.datdir = datdir;
idx.nevdir = nevdir;
idx.exptstarts = sampleons;
idx.exptends = sampleoffs;
idx.exptids = exptids;
   if compareidxs || reindex == 2
       cmpdiffs = CompareGridIdx(idx, oldidx.idx);
       idx.exptdiffs = length(cmpdiffs);
   elseif newf > 0 
       save(idxfile,'idx');
   end
if length(plotexpts)
    GetFigure(checktag);
    hold off;
    for j = 1:length(plotexpts)
        PlotTrialMarks(idx, Expts, plotexpts(j));
        hold on;
    end
elseif plottype == 3
    PlotIdx(idx);
elseif plottype == 1
    PlotExptFiles(idx, Expts, Trials);
end

function PlotExptFiles(idx, Expts, Trials)

for j = 1:length(Expts)
    t(1) = Expts{j}.Header.CreationDate + Expts{j}.Trials(1).Start(1)./(10000 * 60 * 60 * 24);
    t(2) = Expts{j}.Header.CreationDate + Expts{j}.Trials(end).End(end)./(10000 * 60 * 60 * 24);
    plot(t,[1.1 1.1],'r-','linewidth',2);
    text(mean(t), 1.1, Expts{j}.Header.expname,'rotation',90);
    hold on;
end
for j = 1:length(idx.starttime)
    plot([idx.starttime(j) idx.stoptime(j)],[1.9 1.9],'b-','linewidth',2);
    if isnan(idx.stoptime(j))
        text(idx.starttime(j), 1.9, idx.names{j},'rotation',-90);
    else
        text(mean([idx.starttime(j) idx.stoptime(j)]), 1.9, idx.names{j},'rotation',-90);
    end
end

if isfield(Trials,'ExptList')
    id = find([Trials.ExptList.result] == 19); %Canceled expts
    for j = 1:length(id)
        t(1) = Trials.Header.CreationDate + Trials.ExptList(id(j)).start./(10000 * 60 * 60 * 24);
        t(2) = Trials.Header.CreationDate + Trials.ExptList(id(j)).end./(10000 * 60 * 60 * 24);
        plot(t,[1.1 1.1],'k-','linewidth',2);
    end
end
 datetick('x','HH:MM');
 set(gca,'ylim',[1 2]);
   
   
function PlotTrialMarks(idx, Expts,  exptno)
   
        plotsummary = 2;
   id = find(idx.expt == exptno);
   bid = id;
   for j = 1:length(id)
       nfiles{j} = [idx.datdir '/' idx.names{id(j)}];
   end
   if plotsummary == 2 %plot trial start/end
       b = exptno;
       T = Expts{b}.Trials;
       ebstimes = Expts{b}.bstimes;
    for j = 1:length(ebstimes)
%        plot([ebstimes(j) ebstimes(j) eestimes{b}(j) eestimes{b}(j)],[0 1 1 0],'r-');
        plot([ebstimes(j) ebstimes(j)],[0 1],'r-');
        hold on;
    end
    for j = 1:length(T)
        plot([T(j).Start(1) T(j).Start(1)],[-1 -2],'g-');
    end
    for j = 1:length(id)
        bs = idx.bstimes{id(j)} + idx.toff(id(j));
        for k = 1:length(bs)
            plot([bs(k) bs(k)],[0 -1 ],'-');
        end
        text(mean(bs), -1,idx.names{id(j)}(8:10));
        for k = 1:length(idx.evtimes{id(j)})
            t = idx.evtimes{id(j)}(k) + idx.toff(id(j));
            if bitand(1,idx.digin{id(j)}(k))
                plot([t t],[-0.5 -1.5],'g');
            else
                plot([t t],[-0.5 -1.5],'r');
            end
        end
    end
%    aid = find(idx.starttime > Expts{b}.Headerstarts(b) & idx.starttime < ends(b));
    aid = find(idx.toff > Expts{b}.Header.Start & idx.toff < Expts{b}.Header.End);
    for j = 1:length(aid)
        t = (idx.starttime(aid(j))-Expts{b}.Header.CreationDate) * (60 * 60 * 24 * 10000);
        plot([t t],[-2 -3],'r');
        dy = mod(j,5)/5;
        text(t, -2 - dy,idx.names{aid(j)}(8:10));
    end

    for j = 1:length(aid)
        t = idx.toff(aid(j));
        plot([t t],[-2 -3],'b');
        plot(t+idx.ddelay(aid(j)),-2.5,'x');
        dy = mod(j,5)/5;
        text(t, -2 - dy,idx.names{aid(j)}(8:10));
    end
    E = Expts{b};
    requesttime = [];
    if isfield(E,'DigMark')
        for j = 1:length(E.DigMark.times)
            t = E.DigMark.times(j).*10000;
            if bitand(1,E.DigMark.codes(j))
                plot([t t],[-4 -5],'r');
                requesttime = [requesttime t];
            else
                plot([t t],[-4 -5],'g');
            end
        end
    end
    t = E.gridstoreon(1);
    plot([t t],[-3 -4],'r');
    for j = 1:length(E.gridstoreoff)
        id = find(E.gridstoreon < E.gridstoreoff(j));
        if length(id)
        t = E.gridstoreon(id(end));
        plot([t t],[-3 -4],'r');
        end
    end
   end
   
   
   

                

function idx = BuildMatFiles(datdir)

idx = [];
reindex = 0;
nf = 0;

d = dir(datdir);
filenames = {d.name};
newf = 0;
for j = 1:length(d)
    if strfind(d(j).name,'.nev') %this has starttime and Dig Events
        nevfile = [datdir '/' d(j).name];
        matfile = strrep(d(j).name,'.nev','.mat');
        mid = strmatch(matfile,filenames);
        if length(mid)
            agediff = d(j).datenum > d(mid).datenum;
        else
            agediff = 1;
        end
        if reindex || isempty(idx) || isempty(strmatch(d(j).name,idx.names)) || ...
                agediff > 0
            nf = nf+1;
            if  agediff > 0
                nev = openNEV('read','nomat','noparse','nowarning',nevfile);
            else
                nev = openNEV('read','noparse','nowarning','nowaves', nevfile);
            end
            idx.names{nf} = d(j).name;
            idx.expt(nf) = nf;
            if ~isempty(nev.Data.Spikes.Electrode)
            idx.nprobes(nf) = max(nev.Data.Spikes.Electrode);
            end
            idx.toff(nf) = 0;
        end
    end
end
idx.datdir = datdir;
idx.nevdir = datdir;

function PlotIdx(idx)

for j = 1:length(idx.exptstarts)
    plot([idx.exptstarts(j) idx.exptends(j)],[1 1],'-','linewidth',3);
    hold on;
    plot([idx.exptstarts(j) idx.exptstarts(j)],[0 3],'k-');
    text(idx.exptstarts(j),1.2,num2str(j));
end
for j = 1:length(idx.starttime)
    plot([idx.starttime(j) idx.stoptime(j)],[2 2],'r-','linewidth',3);
    hold on;
    text(idx.starttime(j),2.2,num2str(j));
end
set(gca,'ylim',[0 3]);
datetick('x','HH:SS');





function Expts = GetExpts(name, varargin)

args = varargin;
[Trials, Expts, All] = APlaySpkFile(name, 'nospikes',  args{:});
if Trials.Header.Spike2Version < 1.27
    sonid = find(All.Events.codes(:,1) ==48); %'0' = storage on
    soffid = find(All.Events.codes(:,1) ==49); %'1' = storage off
else
    sonid = find(All.Events.codes(:,1) ==49); %'1' = storage on
    soffid = find(All.Events.codes(:,1) ==48); %'0' = storage off
end
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
        Expts{nexp}.gridstoreoff = All.Events.times(soffid(id))./10000;
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


