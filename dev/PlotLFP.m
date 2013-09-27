function result = PlotLFP(Expt,varargin)

lfprange = [30 60 200 10 4];
duration = [];
extraexp{1} = 'st] == 0';
extra.label{1} = 'Blank';
nextra = 1;
forcexvs = 0;
addn=0;
lfpch = 7;
addl=0;
logx = 0;
fillall = 0;
btype = [];
TIMEPLOT = 1;
PWRSPEC = 2;
plottype = TIMEPLOT;

type = Expt.Stimvals.et;
if isempty(strmatch(Expt.Stimvals.e2,'e0'))
    btype = Expt.Stimvals.e2;
end

colors = mycolors;
j = 1;
nvar = nargin - 2;
nmin = 1;
 while j  <= length(varargin)
    str = varargin{j};
    if strncmpi(str,'Type2',5)
        j = j+1;
        if isempty(strmatch(varargin{j},{'e0' type}))
            btype = varargin{j};
        end
    elseif strncmpi(str,'Chan',4)
        j = j+1;
        lfpch = varargin{j};
    elseif strncmpi(str,'power',5)
        plottype = PWRSPEC;
    end
    j = j+1;
 end
 
[Expt, ck] = FixLFPlen(Expt);
 if btype
     if isfield(Expt.Trials,btype)
         bvals = eval(['sort(unique([Expt.Trials.' btype ']));']);
     else
         btype = [];
     end
 end
      plotlfp = 1;
      len = Expt.Header.lfplen;
      if isfield(Expt.Header,'LFPsamplerate')
          ftfrq = (0:len-1)/(len * Expt.Header.LFPsamplerate .* 10);
      else
          ftfrq = (0:len-1)/(len/(10000 * Expt.Header.CRsamplerate));
      end
      result.lfpfrq = ftfrq;
      subplot(2,1,2);
      hold off;
 
      
    fidx = find(ftfrq < lfprange(3));
    fidx = fidx(2:end);
    gidx = find(ftfrq < lfprange(2) & ftfrq > lfprange(1));
    agidx = find(ftfrq < lfprange(1) & ftfrq > lfprange(4));
    bgidx = find(ftfrq < lfprange(4) & ftfrq > lfprange(5));
    cgidx = find(ftfrq < lfprange(2));
    
Expt = FillExpt(Expt, 'st');
if ~isfield(Expt.Trials,type)
    Expt = FillExpt(Expt,type); % just in case
end
stimtype = GetEval(Expt,'st','mode');
[eye, eyes] = GetEval(Expt,'me','mode');

if isfield(Expt.Trials,'ce') 
    if isfield(Expt.Trials,'st')
        sid = find([Expt.Trials.st] > 0);
        isce = sum([Expt.Trials(sid).ce] ~= 0)./length(sid);
        isce = round(isce *4); %at least 25% of stim are correlated, else uncorr is not an interleave
    else
        isce = 1;
    end
    if isce
    nextra = nextra+1;
    extraexp{nextra} = 'ce] == 0';
    extra.label{nextra} = 'Uncorr';
    if (isfield(Expt.Stimvals,'e3') &  strcmp(Expt.Stimvals.e3,'ce')) | ~isempty(strfind(Expt.Header.Name,'ACNC')) | strcmp(Expt.Stimvals.e2,'c1')
        splitextra(nextra) = 1;
    else
        splitextra(nextra) = 0;
    end
    end
end


if isempty(duration)
    duration = min(ck.durs);
end
nch = max(ck.nch);
xvs = sort(unique([Expt.Trials.(type)]));
xvs = xvs(find(~isnan(xvs)));

freq = GetEval(Expt,'tf')/10000;
if(freq > 0)
    period = 1/freq;
else
    period = duration;
end
extraidx = [];

extraidx = find([Expt.Trials.Trial] < 0);
result.excluded = extraidx;
for ie = 1:length(extraexp)
    extra.n(ie) = 0;
    if ~isempty(extraexp{ie})
        %first identify trials that meet condition
        idx = eval(['find([Expt.Trials.' extraexp{ie} ');']);
        %then exclude any that have already been included in earlier condition
        %(makes sure a blank stim is not also counted as uncorr, for example)
        idx = setdiff(idx,extraidx);
        if length(idx) >= nmin & length(idx) < length([Expt.Trials]) * 0.9
            
            extra.x(ie) = mean(xvs);
            extra.n(ie) = length(idx);
            extra.id{ie} = idx;
            extraidx = [extraidx idx];

            fts = abs([Expt.Trials(idx).FTlfp]);
            if nch > 1
                   extra.lfp{ie} = mean(cat(3,Expt.Trials(idx).LFP),3);
                   extra.lfppower{ie} = mean(abs(cat(3,Expt.Trials(idx).FTlfp)),3);
            else
                   extra.lfp{ie} = mean([Expt.Trials(idx).LFP],2);
            end

            
            ft = mean(fts,2);
            aft = mean(fts,1);
            blanklfpwr = mean(ft(gidx));
            ablanklfpwr = mean(ft(agidx));
            bblanklfpwr = mean(ft(bgidx));
            blanklfpse = std(mean(fts(gidx,:),1))/sqrt(length(idx));
            ablanklfpse = std(mean(fts(agidx,:),1))/sqrt(length(idx));
            bblanklfpse = std(mean(fts(bgidx,:),1))/sqrt(length(idx));
            %Calculate mean resp of all non-blank stimuli.
            aidx = setdiff([1:length(Expt.Trials)],extra.id{1});
            if length(aidx)
            oft = mean(abs([Expt.Trials(aidx).FTlfp]),2);
            respf = smooth(oft(fidx)-ft(fidx),3);
            result.lfpower = smooth(oft,5);
            else
                oft = ft;
                respf = smooth(ft(fidx),3);
                result.lfpower = smooth(ft,5);
            end
            [a, peakf] = max(respf);
            j = peakf;
            while(respf(j) > a/5 & j < length(ftfrq) & j < length(respf))
                j = j+1;
            end
            if j > 1 %in case peak is at 1
            maxf = j-1;
            end
            j = peakf;
            while(respf(j) > a/5 & j > 1)
                j = j-1;
            end
            minf = j+1;
            cgidx = minf:maxf;
            cblanklfpwr = mean(ft(cgidx));
            lfpautof = [ftfrq(minf) ftfrq(maxf)];
        else
            blanklfpwr = 0;
            ablanklfpwr = 0;
            bblanklfpwr = 0;
            cblanklfpwr = 0;
            lfpautof = [];
            oft = mean(abs([Expt.Trials.FTlfp]),2);
            result.lfpower = smooth(oft,5);
        end
    end
end
result.extra = extra;
result.Header = Expt.Header;
idx = setdiff(1:length(Expt.Trials),extraidx);

if ~forcexvs
    xvs = sort(unique([Expt.Trials(idx).(type)]));
    xvs = xvs(find(~isnan(xvs)));
    [xvs, mindiff] = RemoveDuplicates(xvs);
end

if exist('ctype','var')
  cvals = eval(['sort(unique([Expt.Trials(idx).' ctype ']));']);
  if strmatch(ctype,'ce')
      cvals = fliplr(cvals);
  end
else
  cvals = NaN;
end


if isempty(btype)
  bvals = NaN;
else
  bvals = eval(['sort(unique([Expt.Trials(idx).' btype ']));']);
end

nlfp = 1;
nfrq = ceil((length(ftfrq)+1)/2);
for nc = 1:length(cvals)
    if ~isnan(cvals(nc))
        cidx  = eval(['find([Expt.Trials.' ctype '] == cvals(nc));']);
    else
        cidx = 1:length(Expt.Trials);
    end
    for ie = 1:length(bvals);
        ix = 1;
        if ~isnan(bvals(ie))
            bidx  = eval(['find([Expt.Trials.' btype '] == bvals(ie));']);
        else
            bidx = 1:length(Expt.Trials);
        end
        step = (max(xvs)-min(xvs))/50;

        % combined files might have .st in one and not antoher
        if isfield(Expt.Trials,'st') & length([Expt.Trials.(type)]) ~= length([Expt.Trials.st])
            for j = 1:length(Expt.Trials)
                if isempty(Expt.Trials(j).st)
                    Expt.Trials(j).st = Expt.Stimvals.st;
                end
            end
        end
  result.linevals(ie) = SetLineval(bvals,bvals(ie));
        for ix = 1:length(xvs);
            x = xvs(ix);
            idx = find(abs([Expt.Trials.(type)]- x) <= mindiff);
            if Expt.Stimvals.st ~= 0
                idx = setdiff(idx,find([Expt.Trials.st] == 0));
            end

            %     idx = find([Expt.Trials.(type)] == x & [Expt.Trials.st] > 0);
            idx = intersect(idx,bidx);
            idx = intersect(idx,cidx);
            idx = setdiff(idx,extraidx);
            if (length(idx) >= nmin | fillall) & (x > 0 | ~logx)
                if x < -999
                    x = SetLineval(xvs,x);
                end
                result.x(ix,ie) = x;
                result.y(ix,ie) = result.linevals(ie);
                result.n(ix,ie) = length((idx));
                ftlfp = cat(3,Expt.Trials(idx).FTlfp);
                ftlfp = squeeze(ftlfp(:,lfpch,:));
                fts = abs(ftlfp);
                ft = mean(abs(ftlfp),2);
                if exist('btype','var')
                    stimlab = [val2str(result.x(ix,ie),type,stimtype,Expt,[]) val2str(result.y(ix,ie),btype,stimtype,Expt,[])];
                else
                    stimlab = [val2str(result.x(ix,ie),type,stimtype,Expt,[])];
                end
                if isempty(idx)
                    fprintf('No Data for %s\n',stimlab);
                else
                        ts = [1:size([Expt.Trials.LFP],1)] .* Expt.Header.LFPsamplerate;
                        result.lfptimes = ts .* 10000;
                        if nch > 1
                            result.lfp(:,ix,ie,:) = mean(cat(3,Expt.Trials(idx).LFP),3);
                            result.lfppower(:,ix,ie,:) = mean(abs(cat(3,Expt.Trials(idx).FTlfp)),3);
                        else
                            result.lfp(ix,ie,:) = mean([Expt.Trials(idx).LFP],2);
                        end
                    if plotlfp == 2 %time domain
                        lfph(nlfp) = plot(ts,mean([Expt.Trials(idx).LFP],2),'color',colors{ix+addn+addl});
                    else
                        lfph(nlfp) = plot(ftfrq(fidx),smooth(ft(fidx),5),'color',colors{ix+addn+addl});
                    end
                    lfplabels{nlfp} = stimlab;
                    nlfp = nlfp + 1;
                    result.lfpwr(ix,ie) = mean(ft(gidx)) - blanklfpwr;
                    result.alfpwr(ix,ie) = mean(ft(agidx)) - ablanklfpwr;
                    result.blfpwr(ix,ie) = mean(ft(bgidx)) - bblanklfpwr;
                    result.clfpwr(ix,ie) = mean(ft(cgidx)) - cblanklfpwr;
                    result.lfpwrs{ix,ie} = mean(fts(gidx,:),1) - blanklfpwr;
                    result.alfpwrs{ix,ie} = mean(fts(agidx,:),1) - ablanklfpwr;
                    result.blfpwrs{ix,ie} = mean(fts(bgidx,:),1) - bblanklfpwr;
                    result.clfpwrs{ix,ie} = mean(fts(cgidx,:),1) - cblanklfpwr;
                    result.lpfsnr = Expt.Header.LFPsnr;
                    hold on;
                    labels{ix} = sprintf('%.1f',x);
                end
            end
            ix = ix+1;
        end
    end
end

result.lfppower = result.lfppower(1:nfrq,:,:,:);
result.ftfrq = result.lfpfrq(1:nfrq);
result.lfpn = result.n;

if isfield(Expt.Header,'LFPtimes')
    result.lfptimes = Expt.Header.LFPtimes .* 10000;
end
subplot(2,1,1);
hold off;
if nch > 1
    plot(result.lfptimes,squeeze(result.lfp(:,:,1,lfpch))');
    if isfield(result.extra,'lfp')
    for j = 1:length(result.extra.lfp)
        hold on;
        plot(result.lfptimes,squeeze(result.extra.lfp{j}(:,lfpch))','k:');
    end
    end
else
plot(result.lfptimes,squeeze(result.lfp)')
end    

  
function Expt = FillExpt(Expt,field)

if ~isfield(Expt.Trials,field)
  for j = 1:length(Expt.Trials)
      Expt.Trials(j).(field) = Expt.Stimvals.(field);
  end
end

function val =SetLineval(bvals,xval)

if xval < -999 %% special interleave
   vmin = min(bvals(find(bvals > -1000)));
   step = mean(diff(bvals(find(bvals > -1000))));
end

if xval == -1003
    val = vmin - step/2;
elseif xval == -1004
    val = vmin - step;
elseif xval == -1001
    val = vmin - step/4;
elseif xval == -1002
    val = vmin - step/2;
else
    val = xval;
end

function str = MakeTitle(Expt,type,type2,stimtype)

str = splitpath(Expt.Header.Name);
if ismember(stimtype,[7,17,19]) %% RDS/Sine
    str = sprintf('%s ic = %.2f',str,GetEVal(Expt,'ic'));
elseif stimtype == 3
    str = [str sprintf(' Or = %.2f SF = %.2f TF = %.2f',GetEval(Expt,'or'),GetEval(Expt,'sf'),GetEval(Expt,'tf'))];
elseif stimtype == 18
    sf = GetEval(Expt,'sf');
    sfb = GetEval(Expt,'f2');
    str = [str sprintf(' SF = %.2f, +- %.2f ',sf,sfb-sf)];
elseif stimtype == 10
    str = [str sprintf('SF = %.2f,%.2f TF = %.2f',GetEval(Expt,'sf'),GetEval(Expt,'f2'),GetEval(Expt,'tf'))];
elseif stimtype == 21 %image
%    str = [str sprintf('SF = %.2f',GetEval(Expt,'sf')) ' \pm ' sprintf('%.2f TF = %.2f  ',GetEval(Expt,'ob'),GetEval(Expt,'tf'))];
     str = [str sprintf('SF = %.2f TF = %.2f  ',GetEval(Expt,'sf'),GetEval(Expt,'tf'))];
end
str = [str sprintf('%.1fx%.1f',GetEval(Expt,'wi'),GetEval(Expt,'hi'))];

function str = val2str(val, type, stimtype,Expt,nfp)

if nargin < 5 | isempty(nfp)
    nfp = 2;
end

if(strmatch(type,'me'))
  if(val == -1)
    str = 'Left';
  elseif(val == 0)
    str = 'Binoc';
  else
    str = 'Right';
  end
elseif(strmatch(type,'sz'))
  str = sprintf('Size %.2f',val);
elseif(strmatch(type,'sd') & ismember(stimtype,[7 17 18]))
    if val == 0
        str = 'RDS';
    elseif val == 1
        str = 'R:grating L:rds';
    elseif val == -1
        str = 'R:rds L:grating';
    elseif val == 2
        str = 'Grating';
    else
        str = sprintf('%.2f',val);
    end
elseif(strmatch(type,'dq'))
    if(val == -1003)
        str = sprintf('sf%.2f',GetEval(Expt,'sf','mode'));
    elseif val == -1004
        str = sprintf('sf%.2f',GetEval(Expt,'f2','mode'));
    else
        str = sprintf('dq=%.2f',val);
    end
elseif(strmatch(type,'dp'))
    if(abs(val) == 1001)
        str = 'Right';
    elseif abs(val) == 1002
        str = 'Left';
    else
        str = sprintf('dp=%.2f',val);
    end
elseif(strmatch(type,'tf'))
  str = sprintf('tf=%.2f',val);
else
  str = sprintf('%s=%.*f',type,nfp,val);
end

function [LFP, ck] = FixLFPlen(LFP)

rate = LFP.Header.LFPsamplerate .* 10000;
for j = 1:length(LFP.Trials)
    lens(j) = size(LFP.Trials(j).LFP,1);
    nch(j)  = size(LFP.Trials(j).LFP,2);
   lfptime = LFP.Trials(j).lfptime;
% Kludege need to think about righ tway to handle lptime = 0 (= missing
% start? 
   if lfptime == 0
       lfptime = LFP.Trials(j).Start(1) - 100;
   end
    first(j) = (LFP.Trials(j).Start(1) - lfptime)./rate;
    tb(j) = lfptime + lens(j) * rate - LFP.Trials(j).End(end);
    dur(j) = (LFP.Trials(j).End(end) - LFP.Trials(j).Start(1))./rate;
end

ck.lens = lens;
ck.nch = nch;
ck.durs = dur;

if length(unique(lens)) == 1
    return;
end

pre = min(first);
tl = min(lens);
    last = (first-pre)+min(lens);
for j = 1:length(LFP.Trials)
    id = round((1+first(j) - pre):last(j));
    if last(j) > lens(j)
        id = id - round(last(j)-lens(j));
    end
    LFP.Trials(j).LFP = LFP.Trials(j).LFP(id,:);    
    %Kludge. Need to recalc this and FTfreq in header
    LFP.Trials(j).FTlfp = abs(fft(LFP.Trials(j).LFP));
    LFP.Header.lfplen = min(lens);
    LFP.Header.LFPtimes = rate * 10000 * ([1:min(lens)] - pre); 
end

function [xvs, mindiff] = RemoveDuplicates(xvs)

mindiff = range(xvs)./200;
if isempty(xvs)
    return;
end
xvs = sort(xvs);
id = find(diff(xvs) > mindiff);
xvs = xvs([1 1+id]);