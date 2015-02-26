function [kernel, details] = BuildBWPsychRCa(Expt, varargin)
% [kernel, details] = BuildBWPsychRC(Expt, varargin)
% Builds Psych and STA kernels in Fourier domain.
% [kernel, details] = BuildBWPsychRC(Expt, 'sprc') to build the STA
% [kernel, details] = BuildBWPsychRC(Expt) to build psych only
% [kernel, details] = BuildBWPsychRC(Expt, 'sprconly') to build the STA wihtout Psych kernel
%                     BuildBWPsychRC(Expt, 'sprconly','fts',details) 
%                supplies Fourier Transforms and seeds to speed up (for resampling)
%                N.B. fts often removed from details before saving, so this typically is
%                done when details has just been calculated
%
%BuildBWPsychRC(kernel, details) to plot a previously calcualted result
%
%BuildBWPsychRC(details) where details is a cell array, containing the
%kernel in a field '.kernel');
%
%BuildBWPsychRC(Expt,'minseedcount',n) changes mininum numnbe of repeated
% %stim required to build kernel (default 200)
%
%BuildBWPsychRC(Expt,'bw',[bw or]) Uses stimuli with orietnation or and
%bandwidth bw to calculate a kernel. Drops minseedcount to 50.
%
%BuildBWPsychRC(Expt,'allstims') Caclulates kernels for each stimulus, and
%then returns an average weighted by the smaller number of choices.
%BuildBWPsychRC(Expt,'byrw') Caclulates kernels separately for each reward
%size
%'allstims' and 'byrw' execute loops when the argument is read, so should
%be the last argument. 
%BuildBWPsychRC(Expt,'parallel, 'allstims') uses a parfor loop
%
%BuildBWPsychRC(Expt,..., 'trialfraction', [a b]) builds a kernel using
%only part of each trial, where a and b are fractions of the duration.
%i.e. [0 0.5] uses only frames from the first half of each trial

rotate = 0;
parallel = 0;
smoothw = 2;
nbins = 6;
seedoffset = 815; %works for lemM079
seedoffsets = 815;
sprc = 0;
delays = [0 300:50:1000];
trialfraction = [0 1];
permute = 0;
spermute = 0;
ftpwr = [];
fts = [];
framerate = 166.7;
minrw = 0;
setrw = 0;
ftseed = 0;
bw = 130;
signalbw = 60;

or = [];
sf = 2;
wi = 1;
kbw = 130;  %BW of stim to use for kernel
kor = 0;
savemode.save = 0;
savemode.seedseq = 0;
minseedcount = [];
rebuild = 0;
pwrargs = {};

nresample = 0;
frameperiod = 166; %default is 60Hz in case
meanimage = 0;
noplot = 0;
plotarg = {};
pargs = {};
pj = 0;
plottype = 'BWRC';
mksigkernel = 0;
predictall = 0;
TrialRange = [];
recorded=0;
paralell = 0;
filterargs = {};


j = 1;
while j <= length(varargin)
    if isstruct(varargin{j})
        details = varargin{j};
    elseif strncmpi(varargin{j},'allstim',5)
        Expt = GetExpt(Expt);
        details.filename = Expt.filename;
        if ~isfield(Expt,'Trials')
            return;
        end
%first build a list of  bw/or conditions        
        ors = unique([Expt.Trials.or]);
        bws = unique([Expt.Trials.ob]);
        nk = 0;
        if sum(bws >= 130)
            nk = nk+1;
            bw(nk) = 130;
            or(nk) = 0;
        end
        na = j;
        for j = 1:length(bws)
            for k = 1:length(ors)
                id = find([Expt.Trials.or] == ors(k) & [Expt.Trials.ob] == bws(j));
                np = sum([Expt.Trials(id).RespDir] == -1);
                nn = sum([Expt.Trials(id).RespDir] == 1);
                if nn > 1 && np > 1 && bws(j) < 120
                    nk = nk+1;
                    bw(nk) = bws(j);
                    or(nk) = ors(k);
                end
            end
        end
        details.allkernels = [];
        if isempty(minseedcount)
            minseedcount = 50;
        end
%Then build the kernels        
        if parallel
            parfor (j = 1:nk)
                Ks{j} = BuildStimKernel(Expt, bw(j), or(j),minseedcount);
            end
        else
            for j = 1:nk
                Ks{j} = BuildStimKernel(Expt, bw(j), or(j),minseedcount);
                nk = 0;
            end
        end
        nk = 0;
        for j = 1:length(Ks)
            if isfield(Ks{j},'nframes')
                nk = nk+1;
                nf(nk) = min(Ks{j}.nframes);
                if nk == 1
                    sumkernel = Ks{j}.kernel .* nf(nk);
                else
                    sumkernel = sumkernel + Ks{j}.kernel .* nf(nk);
                end
                details.allkernels(nk,:,:) = Ks{j}.kernel;
                details.allpsum(nk,:,:) = Ks{j}.psum;
                details.allnsum(nk,:,:) = Ks{j}.nsum;
                details.stimsum(nk,:,:) = Ks{j}.truesum;
                details.or(nk) = or(j);
                details.bw(nk) = bw(j);
                details.nf(nk,:) = Ks{j}.nframes;
                details.zid{nk} = Ks{j}.zid;
                details.Trials(Ks{j}.zid) = Ks{j}.Trials(Ks{j}.zid);
            end
        end
        kernel = sumkernel./sum(nf);
        details.choiceors = ors;
        details.Stimvals = Expt.Stimvals;


        for j = 1:length(details.Trials)
            good(j) = length(details.Trials(j).RespDir);
        end
        details.Trials = details.Trials(find(good));

        good = zeros(1,length(Expt.Trials));
        PlotKernel(kernel,details,plottype);
        if recorded
            details.recorddate=recorddate;
        end
        details.Header = Expt.Header;
        if savemode.save
            if isfield(details.Trials,'rsum')
                exid = [Expt.Trials.id];
                tid = [Expt.Trials.Trial];
                for j = 1:length(details.Trials)
                    id = find(exid == details.Trials(j).id & tid == details.Trials(j).Trial);
                    if length(id) == 1
                        Expt.Trials(id).rsum= details.Trials(j).rsum;
                        good(id) = 1;
                    end
                end
            end
            if sum(good) > length(good)/0.8 && sum(good) < length(good)
                fprintf('Removing incomplete trials %s\n',sprintf('%d ',find(good ==0)));
                Expt.Trials = Expt.Trials(good);
                details.emptyrsum = find(good ==0);
            end
            outname = [details.filename '.mat'];
            Kernel = rmfields(details,{'fts' 'pim' 'nim' 'sigkernel' 'allmean' 'kb' 'tpwr'});
            Kernel.kernel = kernel;
            Kernel.Trials = rmfields(Expt.Trials,'Seedseq');
            Kernel.Header = Expt.Header;
            Kernel.Header.builddate = now;
            Kernel.Header.args = varargin;

            fprintf('Saving Kernel to %s\n',outname);
            save(outname,'Kernel');
        end
        return;
    elseif strncmpi(varargin{j},'bw',2);
        j = j+1;
        kbw = varargin{j}(1);
        kor = varargin{j}(2);
    elseif strncmpi(varargin{j},'byrw',4)
        Expt = GetExpt(Expt);
        details.Trials = Expt.Trials;
        details.filename = Expt.filename;
        [nrw, rws] = Counts([Expt.Trials.rw]); 
        rws = rws(nrw > 100);
        nk = 1;
        for k = 1:length(rws)
            [kernel,b] = BuildBWPsychRCa(Expt,'rebuild','nosigkernel','minseedcount',1,'setrw',rws(k),varargin{1:j-1});
            if isfield(b,'nframes')
            nf(nk) = min(b.nframes);
            details.pkresp{nk} = b.pkresp;
            details.allkernels(nk,:,:) = kernel;
            sumkernel = kernel .* min(b.nframes);
            details.rws(nk) = rws(k);
            details.Trials(b.zid) = b.Trials(b.zid);
            nk = nk+1;
            end
        end
        kernel = sumkernel./sum(nf);
        ors = unique([Expt.Trials.or]);
        details.nf = nf;
        details.Stimvals = Expt.Stimvals;
        details.choiceors = ors;
        for j = 1:length(details.Trials)
            good(j) = length(details.Trials(j).RespDir);
        end
        details.Trials = details.Trials(find(good));
        
        PlotKernel(kernel,details,plottype);
        if savemode.save
            outname = [details.filename '.mat'];
            Kernel = rmfields(details,{'fts' 'pim' 'nim' 'sigkernel' 'allmean' 'kb' 'tpwr'});
            Kernel.kernel = kernel;
            if savemode.seedseq == 0
                Kernel.Trials = rmfields(Expt.Trials,'Seedseq');
            end
            Kernel.Header = Expt.Header;
            if recorded
                Kernel.recorddate=recorddate;
            end
            Kernel.Header.builddate = now;
            save(outname,'Kernel');
        end
        return;
    elseif strncmpi(varargin{j},'blocks',4);
        j = j+1;
        blocks = varargin{j};
        if isfield(Expt,'Header') && isfield(Expt.Header,'BlockStart')
            TrialRange = [];
            bs = Expt.Header.BlockStart;
            bs(length(bs)+1) = length(Expt.Trials);
            for j = 1:length(blocks)
                TrialRange = [TrialRange bs(blocks(j)):bs(blocks(j)+1)];
            end
          fprintf('Using Trials %d to %d\n',min(TrialRange),max(TrialRange));
          break
        end
    elseif strncmpi(varargin{j},'delay',4);
        j = j+1;
        delays = varargin{j};
    elseif strncmpi(varargin{j},'mean',4);
        meanimage = 1;
    elseif strncmpi(varargin{j},'nojack',4);
        pj = pj+1;
        pargs{pj} = varargin{j};
    elseif strncmpi(varargin{j},'bjack',4);
        pj = pj+1;
        pargs{pj} = varargin{j};
    elseif strncmpi(varargin{j},'nosigkernel',6);
        mksigkernel = 0;
    elseif strncmpi(varargin{j},'noplot',4);
        pj = pj+1;
        pargs{pj} = varargin{j};
        plotarg = 'noplot';
        noplot = 1;
    elseif strncmpi(varargin{j},'or',2);
        j = j+1;
        or = varargin{j};
    elseif strncmpi(varargin{j},'minseedcount',6)
        j = j+1;
        minseedcount = varargin{j};
    elseif strncmpi(varargin{j},'parallel',5)
        parallel = 1;
        filterargs = {filterargs{:} 'parallel'};
    elseif strncmpi(varargin{j},'predictall',10)
        predictall = 1;
    elseif strncmpi(varargin{j},'predict',4)
        j = j+1;
        kernel = varargin{j};
        j = j+1;
        details = varargin{j};
        if isempty(ftpwr)
            ftpwr = abs(details.fts);
        end
        pid = find([Expt.Trials.RespDir] == 1 & [Expt.Trials.ob] == 130);
        nid = find([Expt.Trials.RespDir] == -1 & [Expt.Trials.ob] == 130);
        zid = [pid nid];
        ses = [Expt.Trials(zid).Seedseq];
        choices = [Expt.Trials(zid).RespDir];
        [tpwr, ses, useseed] = CalcTrialPower(Expt, ftpwr, zid,pwrargs{:});

        details = PredictChoice(kernel, tpwr, ses, choices, details, varargin{j+1:end});
        ppid = find(details.pkresp >= details.criterion(1));
        pnid = find(details.pkresp < details.criterion(1));
        details.ccp(4) = CalcCP([Expt.Trials(ppid).count],[Expt.Trials(pnid).count]);
        if isfield(details,'pkaresp') %used jacknife
        ppid = find(details.pkaresp >= details.criterion(2));
        pnid = find(details.pkaresp < details.criterion(2));
        details.ccp(5) = CalcCP([Expt.Trials(ppid).count],[Expt.Trials(pnid).count]);
        end
        return;
    elseif strncmpi(varargin{j},'permute',4)
        permute = 1;
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            permute = varargin{j};
        end
    elseif strncmpi(varargin{j},'rwmin',4)
        j = j+1;
        minrw = varargin{j};
    elseif strncmpi(varargin{j},'setrw',4)
        j = j+1;
        setrw = varargin{j};
    elseif strncmpi(varargin{j},'sigkernel',6);
        mksigkernel = 1;        
    elseif strncmpi(varargin{j},'rebuild',4)
        rebuild = 1;
    elseif strncmpi(varargin{j},'resample',4)
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            nresample= varargin{j};
        end
    elseif strncmpi(varargin{j},'save',4)
        savemode.save = 1;
    elseif strncmpi(varargin{j},'spermute',4)
        spermute = 1;
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            spermute = varargin{j};
        end
    elseif strncmpi(varargin{j},'plot',4)
        j = j+1;
        plottype = varargin{j};
    elseif strncmpi(varargin{j},'rotate',4);
        j = j+1;
        rotate = varargin{j} .* pi/180;
    elseif strncmpi(varargin{j},'sprc',4);
            sprc = 1;
            if strncmpi(varargin{j},'sprconly',6)
                sprc = 2;
            end
    elseif strncmpi(varargin{j},'ftpwr',3);
        j = j+1;
        ftpwr = varargin{j};
    elseif strncmpi(varargin{j},'fts',3);
        j = j+1;
        if isfield(varargin{j},'fts') %details strcut
            fts = varargin{j}.fts;
            ftseed = varargin{j}.seedoffset;
        else
            fts = varargin{j};
        end
        for k = 1:size(fts,1)
            trueft{k} = squeeze((fts(k,:,:)));
            ftpwr{k} = abs(trueft{k});
        end
    elseif strncmpi(varargin{j},'rfenvelope',4)
        filterargs = {filterargs{:} varargin{j}};
        j = j+1;
        filterargs = {filterargs{:} varargin{j}};
    elseif strncmpi(varargin{j},'rweight',6);
        plotarg = {plotarg{:} varargin{j}};
    elseif strncmpi(varargin{j},'smoothw',6);
        j = j+1;
        smoothw = varargin{j};
    elseif strncmpi(varargin{j},'trialfraction',9);
        j = j+1;
        trialfraction = varargin{j};
    elseif strncmpi(varargin{j},'Trials',6);
        j = j+1;
        TrialRange = varargin{j};
    elseif strncmpi(varargin{j},'zscore',5);
        pwrargs = {pwrargs{:} varargin{j}};
    elseif strncmpi(varargin{j},'recorded',3);
        j = j+1;
        recorded=1;
        recorddate = varargin{j};    
    end
    j = j+1;
end

if isempty(minseedcount)
    minseedcount = 200; %need to lower this for recording sessions
end

if iscell(Expt) && isfield(Expt{1},'kernel') && isfield(Expt{1},'predchoice')
    smw = smoothw-1;
    if plottype == 'PKOR'
        for j = 1:length(Expt)-smw
            [a,b] = BuildBWPsychRCa(Expt(j:j+smw));
            pkors(j,:) = b.pkor;
            peaks(j) = GetAngle(b.pkor);
        end
        hold off; 
        imagesc(pkors);
        hold on;
        plot(peaks,1:j,'w');
        kernel = a;
        details = b;
    else
    kernel = zeros(size(Expt{1}.kernel));
    n = 0;
    details.choices = [];
    %details.pkresp = [];
    details.choiceors = [];
    allor = [];
    details.pkresp = [];
    F = gcf;
    for j = 1:length(Expt)
        figure(F);
        PlotKernel(Expt{j}.kernel, Expt{j}, plottype, varargin{1:end},'plotrow',[length(Expt) 2 j]);
        kernel = kernel + Expt{j}.kernel * length(Expt{j}.choices);
        details.choices = [details.choices Expt{j}.choices];
        details.choiceors = [details.choiceors Expt{j}.choiceors];
        details.pkresp = [details.pkresp Expt{j}.pkresp];
        allor = cat(1,allor, Expt{j}.pkor);
        n = n+length(Expt{j}.choices);
    end
    GetFigure('MeanKernel');
%    [details.pkor, a] = AddOriSum(kernel');
    details.pkvar = var(kernel(:));
    details.filename = 'Combined';
    [kernel, details] = PlotKernel(kernel, details, plottype);
    end
    return;
elseif iscellstr(Expt)
    for j = 1:length(Expt)
       [a, res{j}] = BuildBWPsychRCa(Expt{j},varargin);
       res{j}.kernel = a;
    end
%    kernel = BuildBWPsychRCa(res);
     kernel = a;
    details = res;
    return;
elseif ischar(Expt)
    ts = now;
    [a, Ex] = ReadSerialFile(Expt,'readexpt','exptype','ORBW');
    %E
    fprintf('Read Serial Took %.2f\n',mytoc(ts));
    if ~isfield(Ex,'Trials') || ~isfield(Ex.Trials,'RespDir')
        return;
    end
    [kernel, details] = BuildBWPsychRCa(Ex,varargin{:});
    details.filename = Expt;
    details.Header = Ex.Header;
    details.Stimvals = Ex.Stimvals;
    if savemode.save
        outname = [Expt '.mat'];
        Kernel = rmfields(details,{'fts' 'pim' 'nim' 'sigkernel' 'allmean' 'kb' 'tpwr'});
        Kernel.kernel = kernel;
        if recorded
            Kernel.recorddate=recorddate;
        end
        Kernel.Header.builddate = now;
%        Kernel.Trials = rmfield(Kernel.Trials,'Seedseq');
        save(outname,'Kernel');
    end
    return;
elseif isfield(Expt,'Trials') && rebuild
    details.plotr = 0;
elseif isfield(Expt,'predchoice') & length(varargin) > 0 & isnumeric(varargin{1}) %its a result file.
    [kernel, details] = PlotKernel(varargin{1}, Expt, plottype, varargin{2:end});
    return;
elseif isfield(Expt,'allkernels') & isfield(Expt,'kernel')
    [kernel, details] = PlotKernel(Expt.kernel, Expt, plottype, varargin{2:end});
    return;
elseif isfield(Expt,'predchoice') & isfield(Expt,'kernel') %its a result file with kernel
    [kernel, details] = PlotKernel(Expt.kernel, Expt, plottype, varargin{1:end});
    return;
elseif isfield(Expt,'predchoice') & isfield(Expt,'fts') %its a result file with powers
    if isempty(ftpwr)
            for j = 1:size(Expt.fts,1)
            ftpwr{j} = squeeze(abs(Expt.fts(j,:,:)));
            end
    end
elseif isnumeric(Expt)
    PlotKernel(Expt, details, plottype);
    return;
else
    if isfield(Expt,'Trials') && ~isfield(Expt.Trials,'imseed')
        [Expt.Trials.imseed] = deal(815);
    end
    details.plotr = 0;
end

px  = 0.0188; %was size on Berlioz for V1 expts, and on Ravel.
if length(or) > 1 && meanimage %build FFT plots for several oris with signal, to get signs right
    [nr,nc] = Nsubplots(length(or));
    for j = 1:length(or)
        nfts = filterim([sf sf/2], [or(j) kbw], wi, 'seedoffset', seedoffset, 'pix2deg', px,'nseeds', 10, 'getft','noplot',filterargs{:});
        sigkernel = fftshift(squeeze(mean(abs(nfts),1)));
        subplot(nr,nc,j);
        imagesc([2 257],[1 256],fliplr(sigkernel'));
        title(sprintf('%.0f',or(j)));
    end
    return;
end

if ~isfield(Expt.Trials,'ob')
    fprintf('%s Not an ORBW file\n',Expt.Header.name);
    kernel = [];
    details.err = 'Not ORBW';
    return;
end

if minrw > 0
    id = find([Expt.Trials.rw] > minrw);
    Expt.Trials = Expt.Trials(id);
elseif setrw > 0
    id = find([Expt.Trials.rw] == setrw);
    Expt.Trials = Expt.Trials(id);
end

sf = Expt.Stimvals.sf;
wi = Expt.Stimvals.wi;
if isfield(Expt.Stimvals,'impx') & Expt.Stimvals.impx < 1
    px = Expt.Stimvals.impx;
end


if isempty(or)
or = Expt.Stimvals.or;
ors = unique([Expt.Trials.or]);
or = min(ors);
end
if isfield(Expt.Trials,'seedoffset')
    id = find([Expt.Trials.ob] > 120);
    seedoffset = median([Expt.Trials.ob]);
end
for j = 1:length(Expt.Trials)
    lens(j) = length(Expt.Trials(j).Seedseq);
    if isempty(Expt.Trials(j).ob)
        Expt.Trials(j).ob = 0;
    end
end
ns = prctile(lens,90);
gid = lens >= prctile(lens,90);
badid = find(lens < prctile(lens,90));
longid = find(lens > ns);
for j = 1:length(longid)
    Expt.Trials(longid(j)).Seedseq = Expt.Trials(longid(j)).Seedseq(1:ns);
end
if kbw < 120
    pid = find([Expt.Trials.RespDir] == 1 & [Expt.Trials.ob] == kbw & gid & [Expt.Trials.or] == kor);
    nid = find([Expt.Trials.RespDir] == -1 & [Expt.Trials.ob] == kbw &gid & [Expt.Trials.or] == kor);
else
    kor = or;
    pid = find([Expt.Trials.RespDir] == 1 & [Expt.Trials.ob] == kbw & gid);
    nid = find([Expt.Trials.RespDir] == -1 & [Expt.Trials.ob] == kbw &gid);
end
spid = find([Expt.Trials.RespDir] == 1 & [Expt.Trials.ob] <= 75 & gid);
snid = find([Expt.Trials.RespDir] == -1 & [Expt.Trials.ob] <= 75 & gid);
RespOris = [median([Expt.Trials(spid).or]) median([Expt.Trials(snid).or]); or or+90];
OriCounts = Counts([Expt.Trials([pid nid]).or]);

if length(nid) < 2 || length(pid) < 2
    kernel = NaN;
    return;
end
[a,b] = Counts(lens);
if min(lens) == max(lens)
    nf = min(lens);
elseif ~isfield(Expt.Stimvals,'nf') || median(lens) == Expt.Stimvals.nf
    nf = median(lens);
else
    fprintf('Forcing nframes to 200');
    nf = 200;
end
errs = find(ismember([Expt.Trials.RespDir],[-1 1]) & lens < nf);
for j = 1:length(errs)
    Expt.Trials(errs(j)).Seedseq(lens(errs(j))+1:nf) = 1;
end
if length(pid) == 0 | length(nid) == 0
    kernel = NaN;
    return;
end
if length(TrialRange)
   pid = intersect(pid,TrialRange);
   nid = intersect(nid,TrialRange);
end
zid = [pid nid];
badid = intersect(badid,zid);
for j = 1:length(badid)
    fprintf('Trials %d only %d seeds\n',badid(j),lens(badid(j)));
end



%some recording files are padded to 200 frames
% only use correct # - important when splitting trial up in fractions
if isfield(Expt.Trials,'nf')
    nf = 1+median([Expt.Trials.nf]);
    for j = 1:length(Expt.Trials)
        if length(Expt.Trials(j).Seedseq) > nf
            xseed = Expt.Trials(j).Seedseq(nf+1:end);
            Expt.Trials(j).Seedseq = Expt.Trials(j).Seedseq(1:nf);
            if sum(xseed > 1) > 0
                fprintf('Trial %d has Seeds after nf\n',j);
            end
        end
    end
end

if isfield(Expt.Header,'frameperiod')
    frameperiod = Expt.Header.frameperiod;
end
if isfield(Expt.Trials,'imseed')
%Trials with empty imseed make a mess. Remove
    for j = 1:length(zid)
        if isempty(Expt.Trials(zid(j)).imseed)
            good(j) = 0;
            Expt.Trials(zid(j)).imseed = 0;
        else
            good(j) = 1;
        end
    end
    zid = zid(good==1);
    pid = pid(find(ismember(pid,zid)));
    nid = nid(find(ismember(nid,zid)));

    [seedcounts, details.seedoffsets] = Counts([Expt.Trials(zid).imseed]);
    if details.seedoffsets(1) == 0 & length(details.seedoffsets) == 2
        seedoffset = details.seedoffsets(2);
        [Expt.Trials(zid).imseed] = deal(seedoffset);
        seedcounts = sum(seedcounts);
        seedoffsets = seedoffset;
    else
        seedoffset = median([Expt.Trials(zid).imseed]);
        seedoffsets = details.seedoffsets;
    end
else
    details.seedoffsets = 815;
    seedcounts = length(zid);
end

if isempty(ftpwr)
    if mksigkernel %calc mean diff between signal containing images, so show what ideal kernel would be
        if signalbw == 0
        sig = min([Expt.Trials.ob]);
        else
            sig = signalbw;
        end
        nfts = filterim([sf sf/2], [or sig], wi, 'seedoffset', seedoffset, 'pix2deg', px,'nseeds', 1000, 'getft','noplot',filterargs{:});
        pfts = filterim([sf sf/2], [or+90 sig], wi, 'seedoffset', seedoffset, 'pix2deg', px,'nseeds', 1000, 'getft','noplot',filterargs{:});
        %Pk is Respdir == 1 (smallest ori) - RespDir==-1 (largest ori)
        sigkernel = fftshift(squeeze(mean(abs(nfts),1)-mean(abs(pfts),1)));
        sigors = [or or+90];
    end
    
    %restrict to trials with enough reps to be worth calculating images
    if length(seedcounts) > 1
        xid = find(seedcounts < minseedcount);
        xxid = [];
        for j = 1:length(xid)
            xxid = [xxid find([Expt.Trials(zid).imseed] == seedoffsets(xid(j)))];
        end
        nid = setdiff(nid, zid(xxid));
        pid = setdiff(pid, zid(xxid));
    end
    zid = [pid nid];
    
    sid = find(seedcounts >= minseedcount);
    if isempty(sid)  %not enough trials with one image set to be useful
        fprintf('Too Few Trials (%d) to be useful\n',max(seedcounts));
        kernel = [];
        details.err = 'Too Few Trials';
        return;
    end
    alltid = [];
    %if they are all 0 , will be set to 815 above
    for j = 1:length(sid)
        seed = seedoffsets(sid(j));
        fts = filterim([sf sf/2], [kor kbw], wi, 'seedoffset', seed, 'pix2deg', px,'nseeds', 1000, 'getft','noplot',filterargs{:});
        for j = 1:size(fts,1)
            trueft{j} = squeeze((fts(j,:,:)));
            pwr{j} = abs(trueft{j});
        end
        tid = find([Expt.Trials(zid).imseed] == seed);
        alltid = [alltid tid];
        if sprc < 2
            [X, ses(:,tid),useseed, truesum] = CalcTrialPower(Expt, pwr, zid(tid),'trialfraction',trialfraction,pwrargs{:});
            for k = 1:length(tid)
                tpwr{tid(k)} = X{k};
                idlist(tid(k)) = Expt.Trials(zid(tid(k))).id;
            end
        end
        ftseed = seed;
    end
elseif sprc < 2
    pwr = ftpwr;
    alltid = find([Expt.Trials(zid).imseed] == seedoffsets(1));
    [tpwr, ses, useseed, truesum] = CalcTrialPower(Expt, pwr, zid,'trialfraction',trialfraction,pwrargs{:});
else
    pwr = ftpwr;
end

if sprc == 2
    tid = find([Expt.Trials(zid).imseed] == ftseed);
    alltid = tid;
    ses = [Expt.Trials(zid).Seedseq];
    useseed = ftseed;
    truesum = zeros(size(pwr{1}));
    truesum = abs(fts(ses(1)));
    [a,b] = find(ses ==0);
    if isempty(a)
        nf = size(ses,1);
    else
        nf = min(a)-1;
    end
    ses = ses(1:nf,:);
    [a,b] = Counts(ses(:),'minval',1);
    for j = 1:length(b)
        truesum = truesum + abs(fts(b(j),:,:)).*a(j);
    end
    truesum = squeeze(truesum)./sum(a);
    tpwr = {};
    useseed = 1:size(ses,1);
end

%if some trials had imseeds not in list (too few), need to remove these
%from zid, pid, nid
zid = zid(alltid);
pid = pid(find(ismember(pid,zid)));
nid = nid(find(ismember(nid,zid)));
choices = [ones(size(pid)) ones(size(nid)).*-1];

if meanimage
    kernel = fftshift(squeeze(mean(pwr)));
    a = kernel'; %gives vertial blob for vertical ori
    imagesc([2 257],[1 256],fliplr(a));
    aor = AddOriSum(a, plotarg);
    [a,b] = GetAngle(aor);
    title(sprintf('OR= %.0f->%.0f (%.0f), bw = %.1f',or,a,b,kbw));
    details.meanor = aor;
    return;
end

sevals = unique(ses);
sevals = sevals(sevals > 0);

if sprc < 2
%kernel is respdir ==1  - respdir == -1
[kernel, details] = CalcPKernel(Expt,pid, nid, pwr, tpwr);
details.choiceors = RespOris;
details.OriCounts = OriCounts;
details = PredictChoice(kernel, tpwr, ses, choices, details,pargs{:});
end
details.choices = choices;
details.zid = zid;
details.impx = px;
details.seedoffset = seedoffset; %for 0 signal
details.seedoffsets = seedoffsets; %all
details.seedcounts = seedcounts; %all
details.Header = Expt.Header;
details.Stimvals = Expt.Stimvals;
details.trialfraction = trialfraction;
details.useseed = useseed;
details.truesum = truesum;

if isfield(Expt.Header,'Name')
    details.filename = Expt.Header.Name;
else
    details.filename = Expt.Header.name;
end

%store enough in Trials so that can build PSFs with ExptPsych
for j = length(Expt.Trials):-1:1
    for f = {'or' 'ob' 'RespDir' 'se' 'Seedseq' 'imseed' 'id' 'Trial'}
        if isfield(Expt.Trials,f{1})
            Trials(j).(f{1}) = Expt.Trials(j).(f{1});
        end
    end
end

if sprc < 2
for j = 1:length(zid)
    Trials(zid(j)).pkresp = details.pkaresp(j);
    Trials(zid(j)).pkrespjk = details.pkresp(j); %done with jacknife
    Trials(zid(j)).rsum = RadialSum(fftshift(tpwr{j})',100);
end
end
if length(TrialRange)
    Trials = Trials(TrialRange);
end

if predictall
    ors = unique([Trials.or]);
    bws = unique([Trials.ob]);
%    bws = bws(bws < 120);
    for j = 1:length(ors)
    for k = 1:length(bws)
        fprintf('Or %.0f, bw %.0f: ',ors(j),bws(k));
        id = find(ismember([Expt.Trials.RespDir],[-1 1]) & [Expt.Trials.ob] == bws(k) & [Expt.Trials.or] == ors(j));
        if length(id)
        nfts = filterim([sf sf/2], [ors(j) bws(k)], wi, 'seedoffset', seedoffset, 'pix2deg', px,'nseeds', 1000, 'getft','noplot',filterargs{:});
        for m= 1:size(fts,1)
            pwr{m} = squeeze(abs(nfts(m,:,:)));
        end
        schoices = [Expt.Trials(id).RespDir];
        [tpwr, ses] = CalcTrialPower(Expt, pwr, id, 'absolute');
        X = PredictChoice(kernel, tpwr, ses, schoices, details,'nojacknife');
        for t = 1:length(id)
            Trials(id(t)).pkresp = X.pkresp(t);
        end
        end
    end
    end
    [tpwr, ses] = CalcTrialPower(Expt, pwr, zid,pwrargs{:});
end

details.Trials = Trials;
details.tpwr = tpwr;
pup = sum(choices ==-1)./length(choices);
rand('state',1001);
checkpermute = 0; %set to 1 for debugging 
for p = permute:-1:1
    id = randperm(length(zid));
    ichoices = choices(id);
%    ichoices = ((rand(1,131) < pup)*2) - 1
    ipid = zid(ichoices==1);
    inid = zid(ichoices==-1); 
 %need to put tpwr in the same order as ichoices, so that ichoices(tid)
 % is [1 1 1 1 1 1 1 ..... -1 -1 -1 -1 -1] becuase CalcPk needs this
 %is removed in the jacknife, the correct image is removed from the kernel.
    tid = [find(ismember(zid,ipid)) find(ismember(zid,inid))];
%    tid = id;
    [ki, a] = CalcPKernel(Expt,ipid, inid, pwr,tpwr(tid));
   
    if checkpermute == 1
    tic;
    for j = 1:length(ichoices)
        if ichoices(j) == 1
            iid = setdiff(ipid, zid(j));
            [ka, b] = CalcPKernel(Expt,iid, inid, pwr,tpwr(tid));
        else
            iid = setdiff(inid, zid(j));
            [ka, b] = CalcPKernel(Expt,ipid, iid, pwr,tpwr(tid));
        end
        ka = fftshift(b.psum-b.nsum);
%        ka = fftshift(ka);
        resp(j) = sum(tpwr{j}(:) .* ka(:));
    end
    crit = prctile(resp,100 * sum(ichoices==1)/length(choices));
    pc(resp <= crit) = -1;
    pc(resp > crit) = 1;
    score = sum(pc==ichoices)./length(choices);
    
  if 0
    k = j;
    for j = 1:length(inid)
        iid = setdiff(inid, inid(j));
        [ka, b] = CalcPKernel(Expt,ipid, iid, pwr,tpwr(tid));
        ka = fftshift(ka);
        resp(j+k) = sum(tpwr{j+k}(:) .* ka(:));
    end
  end
    dur(1) = toc;
    end
    kerrs(p) = a.err;
    ks(p,:,:) = ki;
    xc = corrcoef(kernel(:),ki(:));
    xcs(p,1)= xc(1,2);
    xcs(p,3)= sum(ichoices == choices);
    tic;
    a = PredictChoice(ki,tpwr,ses, ichoices, a);
    dur(2) = toc;
    xc = corrcoef(details.pkresp,a.pkresp);
    xcs(p,2)= xc(1,2);
    pperf(p) = a.predchoice(1);
    pvar(p) = var((ki(:)));
end
if noplot == 0
oldf = gcf;
end
for p = nresample:-1:1
    E = Expt;
        id = Bresample(1:length(zid));
        E.Trials = Expt.Trials(zid(id));
        izid = 1:length(id);
    ichoices = choices(id);
%    ichoices = ((rand(1,131) < pup)*2) - 1
    ipid = find(ichoices==1);
    inid = find(ichoices==-1); 
 %need to put tpwr in the same order as ichoices, so that ichoices(tid)
 % is [1 1 1 1 1 1 1 ..... -1 -1 -1 -1 -1] becuase CalcPk needs this
 %is removed in the jacknife, the correct image is removed from the kernel.
    tid = [find(ismember(izid,ipid)) find(ismember(izid,inid))];
    tid = [id(ipid) id(inid)];
%    tid = id;
    [ki, a] = CalcPKernel(E,ipid, inid, pwr,tpwr(tid));

    x = AddOriSum(ki','noplot',plotarg{:});
    GetFigure('Resamples');
    if p == nresample
        hold off;
    else
        hold on;
    end
    plot(x);
    [kor(p),c,kradius(p)] = GetAngle(x);
end
if nresample && noplot == 0
    hold off;
    plot(kradius.*cos(kor .*pi/180),kradius.*sin(kor .*pi/180),'o');
    hold on; plot(kradius.*cos(pi +(kor .*pi/180)),kradius.*sin(pi+(kor .*pi/180)),'o');
    figure(oldf);
end
if permute
end

%kernel, details] = CalcPKernel(Expt,pid, nid, pwr, tpwr);

if mksigkernel & exist('sigkernel','var');
    details.sigkernel = sigkernel;
    details.sigor = AddOriSum(sigkernel','noplot');
    details.sigors = sigors;
    [details.Trials.idealresp] = deal(NaN);
    for j = 1:length(zid);
        details.Trials(zid(j)).idealresp = details.sigor * details.Trials(zid(j)).rsum';        
    end
    ipid = find([details.Trials.idealresp] > 0);
    inid = find([details.Trials.idealresp] < 0);
    details.idealpkor = mean(cat(1,details.Trials(ipid).rsum))-mean(cat(1,details.Trials(inid).rsum));
end
if sprc < 2
    details.kb = CalcPKernel(Expt,pid(2:end), nid, pwr, tpwr(2:end));
end
if permute
    details.kernelsd = squeeze(std(ks,[],1));
    details.permpred = pperf;
    details.permpkvar = pvar;
end
if nresample
    details.resampleors = kor;
    details.resampler = kradius;
end

for j = 1:size(ses,2)
alltriggers{j} = [];
end
for j = 1:size(sevals,1) %what size(ses,2) < size(sevals,1)
    [ij,ik] = find(ses == sevals(j));
    for k = 1:length(ik)
        alltriggers{ik(k)}= {alltriggers{ik(k)} ij(k)};
    end
end

if sprc < 2
    [details.pkor, a] = AddOriSum(kernel',plotarg);
    details.pkor = a.ork;
    details.sfk = a.sfk;
end
if sprc
    [secounts, sevals] = Counts(ses,'minval',1);
    sxmean = zeros(size(fts,2),size(fts,3));
    for j = 1:length(sevals);
        sxmean = sxmean + squeeze(fts(sevals(j),:,:)) .* secounts(j);
    end
    details.allsxmean = sxmean./sum(secounts);
    details.delays = delays;
    for di = 1:length(delays)
        delay = delays(di);
        spseeds = [];
        nf = size(ses,1);
        k = 1:length(zid);
        for j = 1:length(zid)
            frames = floor((Expt.Trials(zid(j)).Spikes-delay)./frameperiod);
            frames = frames(find(frames >0 & frames <= nf));
            allframes{j} = frames;
            nframes = setdiff(1:nf,frames);
            spkseeds{di,j} = ses(frames,k(j))';
            nospkseeds{di,j} = ses(nframes,k(j))';
            spseeds = [spseeds ses(frames,k(j))']; %seeds associated with spikes
        end

        spsum = zeros(size(pwr{1}));
        spsxsum = zeros(size(pwr{1}));
        [nspk, sv] = Counts(spseeds,'minval',1);
%spsum is spike-triggered amplitude spectrum
%spxsum is spike-triggered Fourier Transform
        for j = 1:length(sv)
            se = sv(j);
            spsum = spsum + nspk(j) .* pwr{se};
            spsxsum = spsxsum + nspk(j) .* trueft{se};
        end
        spmean = spsum./sum(nspk) - details.truesum;
        details.spmean(di,:,:) = fftshift(spmean);
        details.spvar(di) = var(spmean(:));
        sxmean = spsxsum./sum(nspk) - details.allsxmean;
        details.sxmean(di,:,:) = fftshift(sxmean);
        details.sxvar(di) = var(abs(sxmean(:)));
        details.spkseeds = spkseeds;
        details.nospkseeds = nospkseeds;
        nspks(di,sv) = nspk;
    end
    ns = length(ses(:));
    nf = sum(nspks(1,:));
    for ip = 1:spermute
        id = ceil(rand(nf,1) .* ns);
        [pnspk, psv] = Counts(ses(id));
        spsum = pnspk(1) .* pwr{psv(1)};
        for j = 2:length(psv)
            se = psv(j);
            spsum = spsum + pnspk(j) .* pwr{se};
        end
        pm = spsum./nf - details.allmean;
        permvar(ip) = var(pm(:));
    end
    if spermute
        details.permvar = permvar;
    end
    [v,sxt] = max(details.sxvar);
    [v,t] = max(details.spvar);
    nspk = nspks(t,:);
    spk = fftshift(details.spmean(t,:,:)); %%undo fftshift above
    sxspk = fftshift(details.sxmean(t,:,:)); %%undo fftshift above
    for j = 1:length(sevals)
        id = find(sv == sevals(j));
        sid = find(sevals == j);
        if isempty(id) | isempty(sid)
            ps(sevals(j)) = 0;
        else
            ps(sevals(j)) = nspk(sevals(j))./secounts(j); %p(Spike|seed)
            pss(sevals(j),:) = nspks(:,sevals(j))./secounts(j);
        end
        fpwr(sevals(j)) = sum(pwr{sevals(j)}(:) .* spk(:));
        sxpwr(sevals(j)) = abs(sum(sum(fts(sevals(j),:,:) .* sxspk)));
    end
    sevars = var(pss);
%have <= for <, incase the percentile resp is 0;
    [lo,a] = MeanIm(pwr, find(ps <= prctile(ps,30)));
    [hi, b] = MeanIm(pwr, find(ps > prctile(ps,70)));
    details.Imret = [a b];
    for j = 1:length(delays)
        p = pss(:,j);
        lo = MeanIm(pwr, find(ps <= prctile(ps,30)));
        hi = MeanIm(pwr, find(ps > prctile(ps,70)));
        psmean(j,:,:) = fftshift(squeeze(hi-lo));
    end
    [a,b] = sort(fpwr(sevals));  %sorted by power
    k = 1;
    nf = 0;
%if some seeds are missing, than length(fwr) = largest seed, not # of
%seeds
    bincount = sum(secounts)./nbins;
    for j = 1:length(sevals) 
        porder(sevals(b(j))) = k; 
        nf = nf+ secounts(b(j));
        if nf > bincount
            nf = 0;
            k = k+1;
        end
    end
    if sum(porder == k) < 100;
        porder(porder == k) = k-1;
    end
    [a,b] = sort(ps(sevals));
    k = 1; nf = 0;
    for j = 1:length(sevals) 
        forder(sevals(b(j))) = k;
        nf = nf+ secounts(b(j));
        if nf > bincount
            nf = 0;
            k = k+1;
        end
    end
    if sum(forder == k) < 100;
        forder(forder == k) = k-1;
    end
    [a,b] = sort(sxpwr(sevals));%sorted by product with STA
    k = 1; nf = 0;
    for j = 1:length(sevals) 
        sorder(sevals(b(j))) = k;
        nf = nf+ secounts(b(j));
        if nf > bincount
            nf = 0;
            k = k+1;
        end
    end
    if sum(sorder == k) < 100;
        forder(sorder == k) = k-1;
    end
    pscale = std(ps(sevals))./std(fpwr(sevals));
    pfit = polyfit(fpwr(sevals).*pscale,ps(sevals),2);
%nextspk is perdicetd resp from power law fit to (spk pwr kernel .* image) vs rate 
    netspk(sevals) = polyval(pfit,fpwr(sevals).*pscale);
    if noplot == 0
        subplot(2,2,1);
        hold off;
        plot(fpwr(sevals),ps(sevals),'o');
        hold on;
        [a,b] = sort(fpwr(sevals));
        plot(a,netspk(b),'linewidth',2);
        [x,y] = xysmooth(fpwr(sevals),ps(sevals),100);
        plot(x,y,'r');
        [x,y] = xysmooth(fpwr(sevals),ps(sevals),10);
        plot(x,y,'g');
    end
        latency = 500;
    for j = zid
        Expt.Trials(j).Predcount(1) = sum(netspk(Expt.Trials(j).Seedseq(useseed))); %from Spike Triggered power spec
        Expt.Trials(j).Predcount(2) = sum(ps(Expt.Trials(j).Seedseq(useseed))); %from Spike trigged response to seeds
        Expt.Trials(j).count = sum(Expt.Trials(j).Spikes > latency & Expt.Trials(j).Spikes < Expt.Trials(j).dur + latency);
        Expt.Trials(j).sepower = porder(Expt.Trials(j).Seedseq(useseed)); %binned/ordered by image power .* power kernel 
        Expt.Trials(j).sepower(Expt.Trials(j).sepower ==0) = 1;
        Expt.Trials(j).frpower = forder(Expt.Trials(j).Seedseq(useseed));
        Expt.Trials(j).sxpower = sorder(Expt.Trials(j).Seedseq(useseed));
    end
    cfig = gcf;
    details.nspks = sum([Expt.Trials(zid).count]);
    if isfield(Expt.Header,'frameperiod')
        frameperiod = Expt.Header.frameperiod;
        Expt.Stimvals.fz = 10000/frameperiod;
    else
        frameperiod = 166.6;
    end
    details.rc = PlotRevCorAny(Expt,'exp','sepower','Trials',zid,'filltimes',frameperiod,'nmin',10,'box','sdfw',framerate,'noplot');
    netspk = PlotRC(details.rc,'netspk','sdsmooth',5,'timerange',[200 1000],'noplot');
    details.predcount = cat(2,Expt.Trials(zid).Predcount);
    for j = 1:max(porder)
        id = find(porder == j);
        netspkx(j) = mean(fpwr(id));
    end
    if j > length(netspk) %%shouldn't happen,  but can - see ruf1976.c1
        netspk(j) = 0;
    end
    
    for j = zid
        Expt.Trials(j).Predcount(3) = sum(netspk(Expt.Trials(j).sepower)); %from subspace map to binned power groups
    end
    pc = cat(1,Expt.Trials.Predcount)';
    if noplot == 0
    figure(cfig);
    subplot(2,2,2);
    plot(netspkx,netspk,'ro-');
    hold off;
    plot(pc(1,:),[Expt.Trials(zid).count],'o');
    hold on;
    plot(pc(2,:),[Expt.Trials(zid).count],'go');
    plot(pc(3,:),[Expt.Trials(zid).count],'ro');
    xlabel('Predicted Count');
    ylabel('Observed');
    subplot(2,2,3);
    imagesc(fliplr(squeeze(details.spmean(t,:,:))'));
    end
    [details.spkor, spring] = AddOriSum(squeeze(details.spmean(t,:,:))',plotarg);
    if exist('kernel','var')
        xc = corrcoef(kernel(:), details.spmean(t,:,:));
        details.psxcorr(1) = xc(2,1);
        xc = corrcoef(details.pkor, details.spkor);
        details.psxcorr(2) = xc(2,1);
    else
        details.psxcorr = [NaN NaN];
    end
    [a,b,c] = GetAngle(details.spkor);
    title(sprintf('Amp: VarRatio %.2f Pkxc %.3f,%.3f Or %.0f',details.spvar(t)./details.spvar(1),details.psxcorr(1),details.psxcorr(2),a)); 
    details.cp = CalcCP([Expt.Trials(pid).count],[Expt.Trials(nid).count]);
    details.ccp(1) = CalcCP([Expt.Trials(pid).count]-pc(1,ismember(zid,pid)),...
        [Expt.Trials(nid).count]-pc(1,ismember(zid,nid))); %cp after adjust counts using STA of power
    details.ccp(2) = CalcCP([Expt.Trials(pid).count]-pc(3,ismember(zid,pid)),...
        [Expt.Trials(nid).count]-pc(3,ismember(zid,nid))); %adjusted using subspace map to binned power
    details.ccp(3) = CalcCP([Expt.Trials(pid).count]-pc(2,ismember(zid,pid)),...
        [Expt.Trials(nid).count]-pc(2,ismember(zid,nid)));
    if isfield(details,'fixpid')
    details.ccp(4) = CalcCP([Expt.Trials(details.fixpid).count],[Expt.Trials(details.fixnid).count]);
    ppid = find(details.pkresp >= details.criterion(1));
    pnid = find(details.pkresp < details.criterion(1));
    %ccp(5) is the CP based on predicted choices, jacknife
    details.ccp(5) = CalcCP([Expt.Trials(zid(ppid)).count],[Expt.Trials(zid(pnid)).count]);

    details.ccp(6) = CalcCP(pc(1,ppid),pc(1,pnid)); %pred from p(spike|power)
    details.ccp(7) = CalcCP(pc(3,ppid),pc(3,pnid)); %pred from subspace map to binned powers
    details.ccp(8) = CalcCP(pc(2,ppid),pc(2,pnid)); %from spike triggered resp to seeds. 
    end
 
%ccp(9) is the Cp based on overpredicted choices (not using jacknife)
    if isfield(details,'pkaresp') %used jacknife
        ppid = find(details.pkaresp >= details.criterion(2));
        pnid = find(details.pkaresp < details.criterion(2));
        details.ccp(9) = CalcCP([Expt.Trials(zid(ppid)).count],[Expt.Trials(zid(pnid)).count]);
    end
    if noplot == 0
    subplot(2,2,4);
    end
else
    spring.r = 0;
end
    



if sprc < 2
    if noplot == 0
        if rotate
            [xi, yi] = meshgrid([1:wi],[1:wi]);
            xr = (xi-128) .* cos(pi/4) + (yi-128) .* sin(pi/4);
            yr = (yi-128) .* cos(pi/20) - (xi-128) .* sin(pi/20);
        end
        hold off;
        imagesc([2 257],[1 256],fliplr(kernel'));
    end
    [details.pkor, b] = AddOriSum(kernel',plotarg);
    details.pkor = b.ork;
    if b.r < spring.r & noplot == 0
        imagesc([2 257],[1 256],fliplr(kernel'));
        AddOriSum(kernel','maxr',spring.r);
    end
    [a,b,c] = GetAngle(details.pkor);
    if isfield(details,'cp') & noplot == 0
        title(sprintf('PK:CP %.2f, %.2f, %.2f Or%.0f',details.cp,details.ccp(1),details.ccp(2),a));
    end
    details.pkvar = var(kernel(:));
else
    kernel = [];
end
details.fts = fts;


function [details, b] = BuildStimKernel(Expt, bw, or, minseedcount)
    [a,b] = BuildBWPsychRCa(Expt,'rebuild','nosigkernel','bw',[bw or],'minseedcount',minseedcount);
    if isfield(b,'nframes')
        details.nframes = b.nframes;
        details.kernel = a;
        details = CopyFields(details,b,{'psum' 'nsum' 'nframes' 'zid' 'truesum'});
        details.Trials(b.zid) = b.Trials(b.zid);
        details
    else
        details = [];
    end


function [Expt, details] = GetExpt(Expt)
details = [];
if ischar(Expt)
    name = Expt;
    [a,Expt] = ReadSerialFile(name,'readexpt','exptype','ORBW');
    Expt.filename = name;
elseif isfield(Expt,'Header') && isfield(Expt.Header,'loadname')
    Expt.filename = Expt.Header.loadname;
elseif isfield(Expt,'Header') && isfield(Expt.Header,'name')
    Expt.filename = Expt.Header.name;
end
if ~isfield(Expt,'Trials') || isempty(Expt.Trials)
    return;
end
id = find(abs([Expt.Trials.RespDir]) ==1);
Expt.Trials = Expt.Trials(id);



function [tpwr, ses, useseed, isum] = CalcTrialPower(Expt, pwr, zid, varargin)

trialfraction = [0 1]; %use whole trial by default
scoretype = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'absolute',4)
        scoretype = 1;
    elseif strncmpi(varargin{j},'zscore',4)
        scoretype = 2;
    elseif strncmpi(varargin{j},'trialfraction',4)
        j = j+1;
        trialfraction = varargin{j};
    end
    j = j+1;
end

if isempty(zid)
    fprintf('No Trials \n');
end
ses = [Expt.Trials(zid).Seedseq];
nse = find(min(ses') > 0);  %remove trailing zeros
ses = ses(nse,:);
for j = 1:length(Expt.Trials)
    if length(Expt.Trials(j).Seedseq) > max(nse)
    Expt.Trials(j).Seedseq = Expt.Trials(j).Seedseq(1:max(nse));
    end
end


first = 1 + floor(trialfraction(1) * size(ses,1));
last = floor(size(ses,1) * trialfraction(2));
useseed = first:last;
isum = zeros(size(pwr{1}));
for j = 1:length(zid)
    tpwr{j} = MeanIm(pwr,ses(useseed,j));
    isum = isum+tpwr{j};
end
isum = isum./length(zid);
if scoretype == 0 %%express power relative to the mean of all trials.
    for j = 1:length(zid)
        tpwr{j} = tpwr{j} - isum;
    end
elseif scoretype == 2 %%express power as z-score
    tvar = std(cat(3,tpwr{:}),1,3);
    for j = 1:length(zid)
        tpwr{j} = (tpwr{j} - isum)./tvar;
    end
end
ses = ses(useseed,:);


function [mim, good] = MeanIm(pwr, id)

good = 0;
if isempty(id)
    mim = zeros(size(pwr{1}));
    fprintf('Empty list for Meanim');
    return;
end
if length(id) > 5 % bizarerely this is much faster, even if n = 200
mim = pwr{id(1)};
    for j = 2:length(id)
        mim =mim + pwr{id(j)};
    end
good = j;
mim = mim./j;
else
    mim = mean(cat(3,pwr{id}),3);
end


function [kernel, details] = CalcPKernel(Expt, pid, nid, pwr, tpwr)
%new method uses only the mean power spectrum for each Trial 
%makes it easier to deal with combined Expts using different image sets



psum = zeros(size(pwr{1}));
nsum = zeros(size(pwr{1}));
asum = zeros(size(pwr{1}));
pz = 0;
nz = 0;
for j = 1:length(pid)
    if isempty(tpwr{j})
        pz = pz+1;
    else
        psum = psum + tpwr{j}; %in order as zid = [pid nid]
        asum = asum + tpwr{j}; %in order as zid = [pid nid]
    end
end
for j = length(pid)+1:length(pid)+length(nid)
    if isempty(tpwr{j})
        nz = nz+1;
    else
        nsum = nsum + tpwr{j}; %in order as zid = [pid nid]
        asum = asum + tpwr{j}; %in order as zid = [pid nid]
    end
end
pn = length(pid)-pz;
nn = length(nid)-nz;
pim = psum./pn;
nim = nsum./nn;
details.pim = pim;
details.nim = nim;
details.psum = psum;
details.nsum = nsum;
details.allmean = (nsum+psum)./(pn+nn);
kernel = fftshift(pim-nim); %kerel is + choice - -choice
details.nframes = [pn nn];
details.err = NaN;

function [kernel, details] = oldCalcPKernel(Expt, pid, nid, pwr, tpwr)

checktrialpwr = 1;

[a,b] = Counts([Expt.Trials(pid).Seedseq]);
[c,d] = Counts([Expt.Trials(nid).Seedseq]);

psum = zeros(size(pwr{1}));
nsum = zeros(size(pwr{1}));
for j = 1:length(a)
    se = b(j);
    psum = psum + pwr{se} .* a(j);
end
for j = 1:length(c)
    nsum = nsum + pwr{d(j)} .* c(j);
end
pim = psum./sum(a);
nim = nsum./sum(c);
details.pim = pim;
details.nim = nim;
details.psum = psum;
details.nsum = nsum;
details.nframes = [sum(a) sum(c)];

kernel = fftshift(pim-nim); %kerel is + choice - -choice

aim = (nsum+psum)./(sum(a)+sum(c));
details.allmean = aim;

if checktrialpwr
psum = zeros(size(pwr{1}));
nsum = zeros(size(pwr{1}));
for j = 1:length(pid)
    psum = psum + tpwr{j}; %in order as zid = [pid nid]
end
for j = length(pid)+1:length(pid)+length(nid)
    nsum = nsum + tpwr{j}; %in order as zid = [pid nid]
end
pn = length(pid);
nn = length(nid);
tk = psum./pn-nsum./nn;
diffs = fftshift(tk)-kernel;
details.err = sum(abs(diffs(:)));
else
    details.err = NaN;
end


function [kernel, details] = PlotKernel(kernel, details, plottype, varargin)


j = 1; 
plotrow = [];
plotargs = {};
while j <= length(varargin)
    if strncmpi(varargin{j},'plotrow',6)
        j = j+1;
        plotrow = varargin{j};
    elseif strncmpi(varargin{j},'var',3)
    else
        plotargs = {plotargs{:} varargin{j}};
    end
    j = j+1;
end


if ~isfield(details,'pkvar')
    details.pkvar = var(kernel(:));
end
if isfield(details,'pkresp')
    ntrials = length(details.pkresp);
else
    ntrials = 0;
end

if strncmpi(plottype,'RCVar',4)
    subplot(2,2,1);
    hold off;
    [a,t] = max(details.spvar);
    sk = squeeze(details.spmean(t,:,:))';
    imagesc([2 257],[1 256],fliplr(sk));
    [a,sb] = AddOriSum(sk);
    [a,b,c] = GetAngle(a);
    title(sprintf('PWR %d trials. Or %.0f, peak %.0f r%.2f',ntrials,a,b,c));
    subplot(2,2,2);
    hold off;
    plot(details.delays./10,details.spvar./details.spvar(1));
    hold on;
    plot(details.delays./10,details.sxvar./details.sxvar(1),'r');
    legend('pwr','sta');
    subplot(2,2,3);
    hold off;
    [a,t] = max(details.sxvar);
    sk = squeeze(abs(details.sxmean(t,:,:)))';
    imagesc([2 257],[1 256],fliplr(sk));
    [a,sb] = AddOriSum(sk);
    [a,b,c] = GetAngle(a);
    title(sprintf('STA. Or %.0f, peak %.0f r%.2f',a,b,c));
    subplot(2,2,4);
    hold off;
    imagesc(ifft2(squeeze(details.sxmean(t,:,:))));
    
elseif strncmpi(plottype,'BWRCall',4)
    if isfield(details,'spvar') %has spike data
        nr = 2;
        nc = 3;
    else
        nr = 1;
        nc = 2;
    end
    if length(plotrow) > 2
        nr = plotrow(1);
        nc = plotrow(2);
        np = (plotrow(3)-1) * nc;
    else
        np = 0;
    end

    subplot(nr,nc,np+1);
    hold off;

    
    if isfield(details,'spvar') %has spike data
        [a,t] = max(details.spvar);
        sk = squeeze(details.spmean(t,:,:))';
        imagesc([2 257],[1 256],fliplr(sk));
        [a,sb] = AddOriSum(sk);
        [a,b,c] = GetAngle(a);
        title(sprintf('PWR %d trials. Or %.0f, peak %.0f r%.2f',ntrials,a,b,c));
        subplot(nr,nc,4);
        hold off;
        plot(details.delays./10,details.spvar./details.spvar(1));
        hold on;
        plot(details.delays./10,details.sxvar./details.sxvar(1),'r');
        legend('pwr','sta');
        xlabel('ms');

    subplot(2,3,5);
    hold off;
    [a,t] = max(details.sxvar);
    sk = squeeze(abs(details.sxmean(t,:,:)))';
    imagesc([2 257],[1 256],fliplr(sk));
    [a,sb] = AddOriSum(sk);
    [a,b,c] = GetAngle(a);
    title(sprintf('STA. Or %.0f, peak %.0f r%.2f',a,b,c));
    
    subplot(2,3,6);
    hold off;
    imagesc(ifft2(squeeze(details.sxmean(t,:,:))));
    end
    subplot(nr,nc,2+np);    
    hold off;
% center of fftshifted image is 129,129. Fliplr will make this 128,129, so
% change x axis so that 129,129 is plotted at center
    imagesc([2 257], [1 256], fliplr(kernel'));
    [pkor,pb] = AddOriSum(kernel',plotargs{:});
    if ~isfield(details,'pkor')
        details.pkor = pkor;
    end
    if isfield(details,'spvar') %has spike data
        [a,t] = max(details.spvar);
        [a,b,c] = GetAngle(details.pkor);
        title(sprintf('Pk: var %.2f, xc %.3f,%.3f Or%.0f',details.pkvar,details.psxcorr(1),details.psxcorr(2),a));
        details.plotr = max([pb.r sb.r]);
    else
        [a,b,c] = GetAngle(pkor);
        title(sprintf('Pk: var %.2f  Or%.0f N%.0f',details.pkvar,a,ntrials));
        details.plotr = pb.r;
        subplot(nr,nc,1+np);
        hold off; 
        plot(pkor);
        hold on;
        if a < 0
            a = a+180;
        end
        plot([a a],get(gca,'ylim'),':');
        a = max(details.choiceors(:));
        plot([a a],get(gca,'ylim'),'r:');
        plot([a+90 a+90],get(gca,'ylim'),'r:');
        [a,b,c] =  fileparts(details.filename);
        title([b c ' ' num2str(unique(details.choiceors(:))')]);
    end
    if isfield(details,'resampleors')
        GetFigure('Resamples');
        subplot(1,2,1);
        hold off;
        kor = details.resampleors;
        kradius = details.resampler;
        plot(kradius.*cos(kor .*pi/180),kradius.*sin(kor .*pi/180),'o');
        hold on; plot(kradius.*cos(pi +(kor .*pi/180)),kradius.*sin(pi+(kor .*pi/180)),'o');
        subplot(1,2,2);
        id = find(kor < 0);
        kor(id) = kor(id)+180;
        hold off;
        hist(kor);

    end
else
    subplot(2,1,1);
    hold off;
    imagesc([2 257],[1 256],fliplr(kernel'));
    [a,pb] = AddOriSum(kernel');
    [a,t] = max(details.spvar);
    subplot(2,1,2);
    hold off;
    [a,t] = max(details.spvar);
    sk = squeeze(details.spmean(t,:,:))';
    imagesc([2 257],[1 256],fliplr(sk));
    [a,sb] = AddOriSum(sk);
    if sb.r < pb.r
        subplot(2,1,2);
        hold off;
        imagesc([2 257],[1 256],fliplr(sk));
        [a,sb] = AddOriSum(sk,'maxr',pb.r);
        details.plotr = pb.r;
    else
        subplot(2,1,1);
        hold off;
        imagesc([2 257],[1 256],fliplr(kernel)');
        [a,pb] = AddOriSum(kernel','maxr',sb.r);
        details.plotr = sb.r;
    end
    subplot(2,1,1);
    [a,b,c] = GetAngle(details.pkor);
    title(sprintf('Pk: var %.2f, xc %.3f,%.3f Or%.0f',details.pkvar,details.psxcorr(1),details.psxcorr(2),a));
    subplot(2,1,2);
    [a,b,c] = GetAngle(details.spkor);
    title(sprintf('VarRatio %.2f, CP %.2f, %.2f, %.2f  Or%.0f,%.0f',details.spvar(t)./details.spvar(1),...
        details.cp,details.ccp(1),details.ccp(2),a,b));
end
    

function [veco, maxo, r] = GetAngle(ors)
%assumes ors is [0:1180] degrees, as in AddOriSum
    [a,maxo] = max(ors(2:end));
    R = ors(2:end) - min(ors);
    %double the angles - its orientation not diretion. 
    sa = sin([pi/90:pi/90:2*pi]);
    ca = cos([pi/90:pi/90:2*pi]);
    r = sum(R .* (ca + i * sa))./sum(R);
    veco = angle(r).* 90/pi; %halve
    r = abs(r);
    

    
function [sor, details]  = AddOriSum(kernel, varargin)
coss = cos([0:pi/180:pi]);
sins = sin([0:pi/180:pi]);
weightbyr = 0;

fixr = 0;
noplot = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'maxr',4)
        j = j+1;
        fixr = varargin{j};
    elseif strncmpi(varargin{j},'noplot',4)
        noplot = 1;
    elseif strncmpi(varargin{j},'rweight',4)
        weightbyr = 1;
    end
    j = j+1;
end

[sor, details] = RadialSum(kernel, fixr);
cc = details.center;
svar = var(details.sir');
[a,b] = max(svar);
id = find(svar(b:end) > a./20);
r = id(1)+b;
xi = cc(1) + [r .* coss; (r-1) .* coss];
yi = cc(2) + [r .* sins; (r-1) .* sins];
uxi = cc(1) - [r .* coss; (r-1) .* coss];
uyi = cc(2) - [r .* sins; (r-1) .* sins];
details.r = r;

if noplot == 0
    cr = caxis;
    sor = (range(cr) .* (sor - min(sor))./range(sor)) + cr(1);
    set(gca,'xlim',[cc(1)-r cc(1)+r],'ylim',[cc(2)-r cc(2)+r]);
    hold on;
    pcolor(xi,yi,fliplr([sor; sor]));
    pcolor(uxi,uyi,fliplr([sor; sor]));
    [a, b] = max(sor);
    [a, c] = min(sor);
    plot([uxi(1,b) xi(1,b)],[yi(1,b) uyi(1,b)],'r');
    plot([uxi(1,c) xi(1,c)],[yi(1,c) uyi(1,c)],'b');
    shading('flat');
    lr = r*0.8;
    text(cc+lr,cc-lr,'45');
    text(cc+lr,cc+lr,'-45');
    text(cc-lr,cc-lr,'135');
end




function details = PredictChoice(kernel, tpwr, seeds, choices, details, varargin)
%
%tpwr is not the mean power spectrum for each trial. Avoids a lot of
%expensive arracy indexing
%jacknife = 1 is a traidional jacknife, but this is succeptible to bias.
%bias = one half kernel has fewer trials, so jacknkife samples have higher
%variance. 
%jacknife =2 removes the current trial, and one other trial of the opposite
%choice, to eliminate the bias problem.
jacknife = 2;
resptype = 1;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'nojacknife',3)
        jacknife = 0;
    elseif strncmpi(varargin{j},'bjacknife',3)
        jacknife = 1;
    elseif strncmpi(varargin{j},'raw',3)
        resptype = 0; %use raw product, not corrcoeff;
    end
    j = j+1;
end

sk = fftshift(kernel);
if sum(details.nframes) > 100 * length(choices)
    nf = size(seeds,1);
else
    nf = 1; %details.nframes records n stimuli, not frames
end
for j = 1:length(tpwr)
    if resptype == 1
        xc = corrcoef(tpwr{j}(:),sk(:));
        resp(j) = xc(1,2);
    else
        resp(j) = sum(tpwr{j}(:) .*sk(:));
    end
end
cpsum = zeros(size(tpwr{1}));
cnsum = zeros(size(tpwr{1}));
cnn = 0;
cpn = 0;
details.pkresp = resp;

aresp = resp; %copy so resp can be used for jacknife
rnd = rand(length(choices),1);
nf = 1;
if jacknife == 2
%   fsums = sum(pwr(seeds,:,:)); Ntrials * 200 x 256 x 256 = too much
%   memory
% it is crucial here that the order of twpr matches the order of choices
% so watch out when permuting
    pn = details.nframes(1) - nf;
    nn = details.nframes(2) -nf;
    for j = 1:length(choices)
        id = find(choices(j:end) == -1 * choices(j));
            if isempty(id)
                
                id = find(choices(1:j) == -1 * choices(j));
            end
            k = id(ceil(rnd(j)*length(id))); %random trial with the other choice.
        if choices(j) == 1
            psum = details.psum - tpwr{j} .*nf;
            nsum = details.nsum - tpwr{k} .* nf;
        else
            nsum = details.nsum - tpwr{j}.*nf;
            psum = details.psum - tpwr{k}.*nf;
        end
        ski = (psum./pn) - (nsum./nn); %kernel excluding this trial
        if resptype ==1
            xc = corrcoef(tpwr{j}(:),ski(:));
            resp(j) = xc(1,2);
        else
        resp(j) = sum(tpwr{j}(:).*ski(:));
        end
        kr(j) = sum(sk(:).*ski(:));
        ksum(j) = sum(ski(:));
    end
elseif jacknife == 1
%   fsums = sum(pwr(seeds,:,:)); Ntrials * 200 x 256 x 256 = too much
%   memory
% it is crucial here that the order of twpr matches the order of choices
% so watch out when permuting
    for j = 1:length(choices)
        if choices(j) == 1
                pn = details.nframes(1) - nf;
                nn = details.nframes(2);
                psum = details.psum - tpwr{j} .*nf;
                nsum = details.nsum;
        else
            pn = details.nframes(1);
            nn = details.nframes(2) -nf;
            nsum = details.nsum - tpwr{j}.*nf;
            psum = details.psum;
        end
        ski = (psum./pn) - (nsum./nn); %kernel excluding this trial
        if resptype ==1
            xc = corrcoef(tpwr{j}(:),ski(:));
            resp(j) = xc(1,2);
        else
            resp(j) = sum(tpwr{j}(:).*ski(:));
        end
        if resp(j) > aresp(j)
            p = details.psum(:) .* tpwr{j}(:);
            n = details.nsum(:) .* tpwr{j}(:);
            jp = psum(:) .* tpwr{j}(:);
            jn = nsum(:) .* tpwr{j}(:);
        else
            p =0;
        end
    end
else
%calculate projection of each frame onto kernel

end
% set a low criterion to get more + choices
p = sum(choices==-1)./length(choices);
crit = prctile(resp,p.*100);
details.criterion(1) = crit;
correct = (resp >= crit & choices ==1) | (resp < crit & choices == -1);
pc(1) = sum(correct)./length(choices); %performance of jacknife
pc(2) = p.*p + (1-p).*(1-p); %performance of random prediction
if jacknife
    crit = prctile(aresp,p.*100);
    details.criterion(2) = crit;
    details.pkaresp = aresp;
    correct = (aresp >= crit & choices ==1) | (aresp < crit & choices == -1);
    pc(3) = sum(correct)./length(choices); %performance of full kernel
end
details.predchoice = pc;
details.pkresp = resp;
%pkresp is prediected responses from jacknife. 
%pkaresp is predicted responses from raw kernel projeciton
details.jacknife = jacknife;
details.resptype = resptype; 
pid = find(choices ==1);
nid = find(choices == -1);
psum = details.psum;
nsum = details.nsum;
ski = sk;
j = 1;
scores(1) = sum(ski(:).*sk(:));
while scores(j) > 0 & j < 100
    if length(pid) > length(nid)
        [a, id] = max(aresp(pid));
        psum = psum - tpwr{pid(id)};
        pid = setdiff(pid,pid(id));
    else
        [a, id] = min(aresp(nid));
        nsum = nsum - tpwr{nid(id)};
        nid = setdiff(nid,nid(id));
    end
    pn = length(pid)*nf;
    nn = length(nid)*nf;
    ski = (psum./pn) - (nsum./nn); %kernel excluding this trial
    j = j+1;
    scores(j) = sum(ski(:).*sk(:));
end
details.fixpid = pid;
details.fixnid = nid;
