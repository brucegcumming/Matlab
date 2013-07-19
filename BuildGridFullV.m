function result = BuildGridFullV(DATA, Expts, e, probes, varargin)

result = [];

forcebuild =0;
starts = DATA.starts;
ends = DATA.ends;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'force',5)
        forcebuild = 1;
    end
    j = j+1;
end

vid  = [];
        sumtrig = [];
        onetrig = [];
        
        id = find(DATA.idx.expt == e);
        tfirst = (Expts{e}.Trials(1).Start(1) - DATA.preperiod)/10000;
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
        result.askprobes = probes;
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
                if DATA.prebuild
                   rawfile = strrep(filename,'.ns5',['rawp' num2str(p) '.mat']);
                    v = load(rawfile);
                    v.Data = double(v.Data);
                    fprintf('Loaded %s Chspk %d (probe %d)\n',rawfile,p,p);
                else
                v = openNSx('read',filename,['e:' num2str(p)]);
                end
                if DATA.smoothw > 0
                    sm = smooth(v.Data,DATA.smoothw);
                    v.Data = v.Data - sm;
                end
                if ~isfield(v.MetaTags,'Resolution')
                    v.MetaTags.Resolution = v.MetaTags.SamplingFreq;
                end
                if isempty(v.MetaTags.Resolution)
                    v.MetaTags.Resolution = 30000;
                end
                if size(v.MetaTags.VoltsPerDigit,1) < p
                    if isempty(v.MetaTags.VoltsPerDigit)
                        v.MetaTags.VoltsPerDigit =  2.5e-07;
                    end
                    Array = [];
                    vscale = v.MetaTags.VoltsPerDigit(1);
                    FullV.elecxy(1) = NaN;
                    FullV.elecxy(2) = NaN;
                else
                    vscale = v.MetaTags.VoltsPerDigit(p);
                    Array = GetArrayConfig(v);
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
            if isempty(v.MetaTags.ElecLabel)
                FullV.NSlabel = '';
            elseif size(v.MetaTags.ElecLabel,1) < p
                FullV.NSlabel = v.MetaTags.ElecLabel(1,:);
            else
                FullV.NSlabel = v.MetaTags.ElecLabel(p,:);
            end
            FullV.samper = 1./30000.237;
            if ~isempty(xid)
               FullV.V = FullV.V(vid);
            elseif DATA.chopfile && length(FullV.blkstart) ==1 && isempty(xid)
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
               FullV.preperiod = DATA.preperiod;
               FullV.postperiod = DATA.postperiod;
               FullV.chopgap = DATA.gaplen;
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
            if DATA.submains == 2
                if isempty(mcycle)
                    [mainsfullv, mcycle] = RemoveMains(FullV, DATA.mainstimes,[],'calconly');
                end
            elseif DATA.submains %if 2, just calculated for first channel
                if isempty(mcycle)
                    [FullV, mcycle] = RemoveMains(FullV, DATA.mainstimes,[]);
                else
                    FullV = RemoveMains(FullV, DATA.mainstimes, mcycle);
                end
            end
            if DATA.buildmean == 3
                if isempty(sumv)
                    sumv = FullV.V./std(FullV.V);
                else
                    sumv(1:length(FullV.V)) = sumv(1:length(FullV.V)) + FullV.V./std(FullV.V);
                end
            end
            FullV.sumscale = 0; %set below if mean is subtracted.    
            xscale = v.MetaTags.Resolution./max(abs(FullV.V));
            FullV.highpass = DATA.highpass;
            if isnan(DATA.highpass)
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
            result.buildtime(p) = FullV.savetime;
            t = getCurrentTask();
            res.workerid = t.ID;
            fprintf('Writing P%d(%d) %s (%s) (Worker %d took %.2f at %s)',p,FullV.chspk,outfile,outvarname,t.ID,mytoc(ts),datestr(now));
            ts = FullV.savetime;
            SaveFullV(outfile, FullV);
            fprintf('+%.2f)\n',mytoc(ts));
        end
        if DATA.buildmean == 3
            if isempty(sumv)
                ts = now;
                sumv = ProcessGridFullV(DATA.idx.datdir,'expts',e,Trials,Expts,DATA.idx,'buildmean','probes',probes);
            end
            sumsq = sumv*sumv';
            for p = probes
                outfile = [DATA.idx.datdir '/Expt' num2str(e) '.p' num2str(p) 'FullV.mat'];
                X = load(outfile);
                FullV = CopyFields(FullV,X.FullV);  %Not sure if any new fields added to FullV after saving.....
                FullV.V = double(X.FullV.V);
                FullV.sumscale = (sumv(1:length(FullV.V))*FullV.V')./sumsq;
                FullV.V = FullV.V - sumv(1:length(FullV.V)) .* FullV.sumscale;
                FullV.V = int16(FullV.V);
                FullV.savetime(2) = now;
                SaveFullV(outfile, FullV);
            end
            save(sprintf('%s/Expt%dFullVmean.mat',DATA.idx.datdir,e),'sumv')
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
    

function took = SaveFullV(name, FullV)

tic;
save(name,'FullV');
took = toc;



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

