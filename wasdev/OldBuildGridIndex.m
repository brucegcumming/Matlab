function idx = BuildGridIndex(name, Expts, varargin)
%idx = BuildGridIndex(name, Expts, ...)
%Build an index of which neV files match Expt .mat files
%name is a filename or direcotry where the .mat file lives (used to
%construct path for where .nev files live
% Expts is a cell array of Expts, as returned by APlaySpkFile
% BuildGridIndex(filename, [], 'reindex') to rescan nev files

    idx = [];
   datdir = 'F:/Utah/jbe/';
   reindex = 0;
   plotexpts = [];
   checktag = 'TrialCheck';
   preperiod = 5000;
   postperiod = 5000;
   maxgap = 10000; %default is to chop if gap > 1 sec
   DigMark = [];
   plottype = 0;
   Trials = [];
   
j = 1;
while j <= length(varargin)
    if isfield(varargin{j},'ExptList')
        Trials = varargin{j};
    elseif strncmpi(varargin{j},'digmark',5)
        j = j+1;
        DigMark = varargin{j};
    elseif strncmpi(varargin{j},'plotfiles',5)
        plottype = 1;
    elseif strncmpi(varargin{j},'plotexpts',5)
        j = j+1;
        plotexpts = varargin{j};
    elseif strncmpi(varargin{j},'plotidx',5)
        plottype = 3;
    elseif strncmpi(varargin{j},'reindex',5)
        reindex = 1;
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
   plotsummary = 2;

   if isempty(strfind(path,'BlackRock'))
       path(path,'/bgc/bgc/matlab/BlackRock');
   end
   idxfile = [datdir '/FileIdx.mat'];
   
   if ischar(Expts)
   end
   
   if exist(idxfile,'file') & ~reindex
       a = load(idxfile);
       idx = a.idx;
       idx.datdir = datdir;
       nidx = length(idx.names);
       if plottype == 1
           PlotExptFiles(idx, Expts, Trials);
       end
       return;
   else
       nidx = 0;
   end
   
   
   if isempty(Expts)
       if isdir(name)
           idx = BuildMatFiles(datdir);
           return;
       else
           Expts = GetExpts(name);
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
       nson = length(sampleon);
       sampleons = [sampleons sampleon(1)];

       id = find(Expts{j}.DigMark.codes ==2);
       
 %Nov 2012 Onwards, only allow 1 ns5 per expt. Much easier to find. 
       if isempty(id)
           fprintf('Expt %d Missing Sample off marker\n',j);
           sampleoff = Expts{j}.Header.CreationDate + Expts{j}.gridstoreoff(1)./(60 *60 *24 *10000);
           offtimes = [offtimes Expts{j}.gridstoreoff(1)]
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
   d = dir([datdir]);
   filenames = {d.name};
   %d = dir([datdir '/*.nev']);
   newf = 0;
   nnev = 0;
   usenev = 1;
   for j = 1:length(d)
       if ~isempty(strfind(d(j).name,'.nev')) && usenev %this has starttime and Dig Events
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
                   nevoff =  nev.Data.SerialDigitalIO.TimeStampSec(id(end));
                   if nevoff < 1
                       fprintf('Nominally short file %.2f\n', nevoff);
                       nevoff =  nev.Data.SerialDigitalIO.TimeStampSec(end);
                       if isfield(nev.Data.Spikes,'TimeStamp') && ~isempty(nev.Data.Spikes.TimeStamp)
                           lastspk = double(nev.Data.Spikes.TimeStamp(end))./nev.MetaTags.TimeRes;
                           nevoff = max([nevoff lastspk]);
                       end
                       t = tstart+nevoff./(60 * 60 * 24);
                       xid = find(sampleoffs) > t;
                       fprintf('Missing ~ %.2f sec of data\n',sampleoffs(id(1))-1);
                   end
                   t = tstart+nevoff./(60 * 60 * 24);
                   tstop = t;
                   [a,c] = min(abs(t - sampleoffs));
                   b = offids(c);
                   tdiff = offtimes(c)-nevoff; %difference in sec
                   cpuclockdiff = (tstop-sampleoffs(c)).*(60*60*24);
                   requesttime = sampleons(c) * 10000;
               end
               if requesttime > 0
                   nbs = length(ebstimes{b});
                   E = Expts{b};
                   idx.nevstarts(nnev) = tstart;
                   nevfiles{nnev} = d(j).name;
               else
                   nbs = 0;
                   ns = 0;
                   needed(j) = 0;
               end
               if ns && nbs
%normally, the first stimon is a trigger send by binoc, and this is not in the ns5 file
%So ebstimes{b} should be one longer than bstimes. If this is true and
%difftime looks sensible,use it
                   if length(ebstimes{b}) == length(bstimes)+1 && ...
                       abs(ebstimes{b}(2)-(bstimes(1)+tdiff*10000)) < 40000
                       bsoff = 2;
                       bsstd = std(ebstimes{b}(2:end)-(bstimes+tdiff*10000));
                       diffs = ebstimes{b} - (bstimes(1)+tdiff*10000);
                       if bsstd > 1000
                           [a,bsoff] = min(abs(diffs));
                       else
% if bssted is low, diffs(bsoff)should always meet criterion below?
%                           diffs(bsoff) = 0;
                       end
                   else
                       diffs = ebstimes{b} - (bstimes(1)+tdiff*10000);
                       [a,bsoff] = min(abs(diffs));
                       if bsoff == 1 & length(ebstimes{b}) > length(bstimes)%suspicious - check that 2 is not better
                           bstd(1) = std(ebstimes{b}(1:ns)-(bstimes+tdiff*10000));
                           bstd(2) = std(ebstimes{b}(2:ns+1)-(bstimes+tdiff*10000));
                           if diff(bstd) < 0
                               bsoff = 2;
                           end
                           bsstd = bstd(bsoff);
                       else
                           if ns <= nbs
                               bsstd = std(ebstimes{b}(1:ns)-(bstimes(1:ns)+tdiff*10000));
                           else
                               bsstd = std(ebstimes{b}(1:nbs)-(bstimes(1:nbs)+tdiff*10000));
                           end
                       end
                   end
               end

               if ns && length(estimes) && nbs && ...
                       (abs(diffs(bsoff)) < 40000 || bsstd < 100) && abs(cpuclockdiff) < 4
                   nidx = nidx+1;
                   newf = newf+1;
                   idx.names{nidx} = d(j).name;

                   idx.expt(nidx) = b;
                   idx.tdiff(nidx) = (tstop-sampleoffs(c)).*(60*60*24);

                   last = min([bsoff+length(bstimes)-1 length(ebstimes{b})]);
                   
%                   extimes = ebstimes{b}(bsoff:last);
%                   xc = corrcoef(diff(bstimes),diff(extimes));
%toff is the timestamp in spike2 associated with t = 0 in the Cerebrus file
                   idx.toff(nidx) = ebstimes{b}(bsoff)-bstimes(1);
                   nt = length(bstimes);
                   if length(ebstimes{b}) >= length(bstimes)+bsoff-1;
                       bsdiffs = ebstimes{b}(bsoff:bsoff+nt-1)-bstimes;
                       idx.bsstd(nidx) = std(bsdiffs);
                       idx.bsdrift(nidx,:) = bsdiffs([1 end]);
                   else
                       idx.bsstd(nidx) = NaN;
                   end
                   [a, subid] = min(abs(idx.toff(nidx)-requesttime));
%delay is best estimate of the delay between the request from spike2 and the file opending
%in Cerebrus. Based where t=0 falls in spike2, relative to requesttime (dig
%marker in spike2). ddelay (below) is time of first recorded event. Can be
%later if all timing pulses are missed
                   idx.delay(nidx) = idx.toff(nidx)-requesttime(subid);
                   idx.firstbs(nidx) = bsoff;
                   idx.start(nidx) = tstart + bstimes(1)./(10000 * 60 * 24); %datenum
                   idx.end(nidx) = tstart + estimes(end)./ (10000 * 60 * 24);
                   idx.nt(nidx) = length(estimes);
                   idx.bstimes{nidx} = bstimes;
                   idx.evtimes{nidx} = nev.Data.SerialDigitalIO.TimeStampSec .* 10000;
                   idx.digin{nidx} = bitand(7,nev.Data.SerialDigitalIO.UnparsedData);
                   idx.digdt(nidx) = diffs(bsoff);
                   cid = find(cuts(:,1) > idx.toff(nidx) & cuts(:,1) < idx.toff(nidx)+nevoff*10000);
                   idx.cuts{nidx} = cuts(cid,:);
                   needed(j) = 1;
               else
                   nidx = nidx+1;
                   idx.expt(nidx) = 0;
                   idx.tdiff(nidx) = tdiff;
                   idx.bsstd(nidx) = 0;
                   if ns
                   idx.digdt(nidx) = diffs(bsoff);
                   else
                   idx.digdt(nidx) = NaN;
                   end
                   needed(j) = 0;
               end
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
               idx.bsoff(nidx) = diffs(bsoff);
               if length(nev.Data.SerialDigitalIO.UnparsedData) 
                   idx.lastdio(nidx) = bitand(7,nev.Data.SerialDigitalIO.UnparsedData(end));
               end
           end
       end
   end
   if isempty(idx)
       fprintf('Missing NeV Data Filesin %s\n',name);
       return;
   end
   
   missing = setdiff(1:length(Expts),idx.expt);
   if length(missing)
       idx.missing = missing;
       warndlg(sprintf('Missing Expts %s',sprintf('%d ',missing)),'Expts Missing');
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
   
   
   if newf
       save(idxfile,'idx');
   end
idx.datdir = datdir;
idx.exptstarts = sampleons;
idx.exptends = sampleoffs;

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


