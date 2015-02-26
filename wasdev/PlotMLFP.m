function [lfpm, result] = PlotMLFP(Expt, varargin)
%[fpm, details] = PlotMLFP(Expt, ..)
%Plot LFP Data for Arrays
%Default is a line plot, one for each probe
%
%fpm = PlotMLFP(Expt, 'stack')
%fpm = PlotMLFP(Expt, 'image')
%fpm = PlotMLFP(Expt, 'split')
%fpm = PlotMLFP(Expt, 'expts', {'jv' 'sf' 'st' 'tf' 'nsfs'},'lines','ftpwr');
% builds pwr spectrum for each combination of expt vars. 
%
% see also PlotLFPpwr

plottype = 1;
nfreq = 500;
checkdim = 0;
minlen = 0;
chscale = 0;
needft = 0;
stackoff = 0.05;
plotlines = 0;
calcresp = 0;
timerange = [];
freqs = 1:nfreq;
yvals = [];
xvals = [];
et = [];
eb = [];
ec = [];
ed = [];
result = [];

if isfield(Expt.Header,'exptvars')  && ~isempty(strfind(Expt.Header.Options,'+exm'))
    extypes = split(Expt.Header.exptvars,',');
    et = extypes{1};
    if length(extypes) > 1
        eb = extypes{2};
    end
else
    extypes = {};
end
probes = [];
styles = {'-' '--' '-.' ':' '-' '--' '-.' ':' '-' '--' '-.' ':' '-' '--' '-.' ':'};
styles = [styles(:); styles(:)];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'calcresp',4)
        calcresp = 1;
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            calcfreq = varargin{j};
        end
    elseif strncmpi(varargin{j},'image',4)
        plottype = 2;
    elseif strncmpi(varargin{j},'check',4)
        checkdim = 1;
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            minlen = varargin{j};
        end
    elseif strncmpi(varargin{j},'csd',3)
        plottype = 7;
    elseif strncmpi(varargin{j},'expt1',5)
        j = j+1;
        et = varargin{j};        
    elseif strncmpi(varargin{j},'expts',5)
        j = j+1;
        extypes = varargin{j};
        et = extypes{1};
    elseif strncmpi(varargin{j},'expt2',5)
        j = j+1;
        extypes{1} = varargin{j};        
    elseif strncmpi(varargin{j},'expt3',5)
        j = j+1;
        extypes{2} = varargin{j};        
    elseif strncmpi(varargin{j},'expt4',5)
        j = j+1;
        extypes{3} = varargin{j};        
    elseif strncmpi(varargin{j},'select',5)
        j = j+1;
        id = eval(['find([Expt.Trials.' varargin{j} ');']);
        Expt.Trials = Expt.Trials(id);
    elseif strncmpi(varargin{j},'scale',4)
        chscale = 20;
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            chscale = varargin{j};
        end
    elseif strncmpi(varargin{j},'lines',4)
        plotlines = 1;
    elseif strncmpi(varargin{j},'split',4)
        plottype = 8;
    elseif strncmpi(varargin{j},'probes',4)
        j = j+1;
        probes = varargin{j};
    elseif strncmpi(varargin{j},'stack',4)
        plottype = 5;
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1; 
            stackoff = varargin{j};
        end
    elseif strncmpi(varargin{j},'ftpwr',4)
        plottype = 9;
        needft = 1;
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1; 
            freqs = varargin{j};
        end
    elseif strncmpi(varargin{j},'ftimage',4)
        plottype = 4;
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1; 
            nfreq = varargin{j};
        end
    elseif strncmpi(varargin{j},'timerange',5)
            j = j+1; 
            timerange = varargin{j};
    elseif strncmpi(varargin{j},'ftstack',4)
        plottype = 6;
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1; 
            nfreq = varargin{j};
        end
                    freqs = 1:nfreq;
    elseif strncmpi(varargin{j},'ft',2)
        plottype = 3;
                if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1; 
            nfreq = varargin{j};
                end
    elseif strncmpi(varargin{j},'xvals',5)
            j = j+1; 
            xvals = varargin{j};
    elseif strncmpi(varargin{j},'yvals',5)
            j = j+1; 
            yvals = varargin{j};
    end
    j = j+1;
end

if Expt.Header.rc
    return;
end

if isfield(Expt.Trials,'lfpstart') %Utah grid
for j = 1:length(Expt.Trials)
    if ~isnan(Expt.Trials(j).lfpstart);
        good(j) = true;
    end
end
Expt.Trials = Expt.Trials(good);
end


if isempty(et)
    et = GetEval(Expt,'et');
end
if isempty(eb)
    eb = GetEval(Expt,'e2');
end
if isempty(ec)
    ec = GetEval(Expt,'e3');
end


if ~isfield(Expt.Trials,et)
    [Expt.Trials.(et)] = deal(0);
end

if isempty(yvals)
    if ~strcmp(eb,'e0') && isfield(Expt.Trials,eb)
        yiv = [Expt.Trials.(eb)];
    else
        yiv = ones(size(Expt.Trials));
    end
    yvals = unique(yiv);
end

if isempty(extypes)
    extypes{1} = et;
    if isfield(Expt.Trials,eb)
        extypes{2} = eb;
    elseif isfield(Expt.Trials,'st')
        extypes{2} = 'st';
    else
        extypes{2} = '';
    end
    if isfield(Expt.Trials,ec)
        extypes{3} = ec;
    else
        extypes{3} = '';
    end
end

for j = 1:length(extypes)
    if ~isempty(extypes{j})
        if ~isfield(Expt.Trials,extypes{j})
            fprintf('No Field %s in Trials\n',extypes{j});
            return;
        end
        if ischar(Expt.Trials(1).(extypes{j}))
            vals{j} = {Expt.Trials.(extypes{j})};
        else
        vals{j} = [Expt.Trials.(extypes{j})];
        end
        uvals{j} = unique(vals{j});
        env(j) = length(uvals{j});
    end
end
    

if ~isnan(et)
    xiv = [Expt.Trials.(et)];
    allx = unique(xiv);
else
    [Expt.Trials.dummy] = deal(1);
    et = 'dummy';
    xiv = [Expt.Trials.dummy];
end

if isempty(xvals)
    xvals = unique(xiv);
else
    xvals = allx(xvals);
end



if checkdim
    nch = 24;
    for j = 1:length(Expt.Trials)
        lens(j) = size(Expt.Trials(j).LFP,1);
        nchs(j) = size(Expt.Trials(j).LFP,2);
    end
    nch = max(nchs);
    if minlen == 0
        len = min(lens);
        good = 1:length(Expt.Trials);
    else
        if ismember(plottype,[3 4 6 9])
            needft = 1;
        end
        good = find(lens > minlen & nchs == nch);
        len = minlen;
    end
    Expt.Trials = Expt.Trials(good);
    for j = 1:length(Expt.Trials)
        Expt.Trials(j).LFP = Expt.Trials(j).LFP(1:len,:);
        if chscale > 1
            Expt.Trials(j).LFP(:,17:nch) = Expt.Trials(j).LFP(1:len,17:nch)./chscale;    
        end
        if needft
            Expt.Trials(j).FTlfp = fft(Expt.Trials(j).LFP);
        end
    end
    xiv = [Expt.Trials.(et)];
    xvals = unique(xiv);
else
    if isfield(Expt.Trials,'LFP')
        len = size(Expt.Trials(1).LFP,1);
    else
        len = 0;
    end
        if ismember(plottype,[3 4 6])
            needft = 1;
        end
end


      if isfield(Expt.Header,'LFPsamplerate')
          ftfrq = (0:len-1)/(len * Expt.Header.LFPsamplerate .* 10);
          if isfield(Expt.Trials,'lfptime');
              presmp = median([Expt.Trials.Start] - [Expt.Trials.lfptime]) ./( 10000 * Expt.Header.LFPsamplerate);
          else
              presmp = mean([Expt.Trials.Start] - [Expt.Trials.ftime]) .* Expt.Header.LFPsamplerate;
          end
          lfpt = ([1:len] - presmp) .* Expt.Header.LFPsamplerate;
      else
          ftfrq = (0:len-1)/(len/(10000 * Expt.Header.CRsamplerate));
          presmp = mean([Expt.Trials.Start] - [Expt.Trials.ftime]) .* Expt.Header.CRsamplerate;
          lfpt = ([1:len] - presmp) .* Expt.Header.CRsamplerate;
      end

LFP = cat(3,Expt.Trials.LFP);
if isfield(Expt.Trials,'FTlfp') & needft
PWR = cat(3,Expt.Trials.FTlfp);
elseif needft
    for j = 1:length(Expt.Trials)
        Expt.Trials(j).FTlfp = abs(fft(Expt.Trials(j).LFP));
        PWR = cat(3,Expt.Trials.FTlfp);
    end
else
    PWR = LFP;
end


for j = 1:length(extypes)
    if ~isempty(extypes{j})
    cx{j} = 1;
    end
end

allids = {};
if length(env) == 1
    env(2) = 1;
end
for ii = 1:prod(env(1:end))
    [cx{:}] = ind2sub(env,ii);
    idx = 1:length(vals{1});
    for k = 1:length(cx)
        if iscellstr(vals{k})
            id = find(strcmp(uvals{k}{cx{k}},vals{k}));
        else
            id = find(vals{k} == uvals{k}(cx{k}));
        end
        idx = intersect(idx,id);
    end
    j = cx{1};
    [j,k] = ind2sub([env(1) prod(env(2:end))],ii);
    n = sum(ismember(idx,cat(2,allids{:})));
    if n == 0
        if isempty(idx)
            idx = [];
        end
        allids{k,j} = idx;
    elseif n < length(idx)
        idx = [];
    else
        idx = [];
    end
    alln(k,j) = length(idx);
    if ~isempty(idx)
        %    for k = 1:size(LFP,2);
        lfpm(k,j,:,:) = squeeze(mean(LFP(:,:,idx),3))';
        lfpft(k, j,:,:) = squeeze(mean(abs(PWR(:,:,idx)),3))';
        %    end
    end
end

result.ids = allids;
if isfield(Expt.Header,'LFPfreq')
    result.freqs = Expt.Header.LFPfreq(freqs);
else
    period = max(Expt.Header.LFPtimes)-min(Expt.Header.LFPtimes);
    result.freqs = [1:length(Expt.Header.LFPtimes)]./period;
    result.freqs = result.freqs(freqs);
end

if isempty(probes)
    probes = 1:size(lfpm,3);
end

colors = mycolors;
if plottype == 1
    for j = 1:size(LFP,2)
        subplot(5,5,j);
        hold off;
        for k = 1:length(xvals)
            plot(squeeze(lfpm(k,j,:)),'color',colors{k});
            hold on;
        end
    end
elseif plottype == 3
    for j = 1:size(LFP,2)
        subplot(5,5,j);
        hold off;
        for k = 1:length(xvals)
            plot(ftfrq(1:nfreq),squeeze(lfpft(k,j,1:nfreq)),'color',colors{k});
            hold on;
        end

    end
elseif plottype == 2 'image'
    [nr,nc] = Nsubplots(length(xvals));
    pmax = max(lfpm(:));
    pmin = min(lfpm(:));
    for k = 1:length(xvals)
        subplot(nr,nc,k);
        imagesc(Expt.Header.LFPtimes,size(lfpm,2),squeeze(lfpm(k,:,:)));
        title(sprintf('%.2f',xvals(k)));
        caxis([pmin pmax]);
    end
    if ~isempty(timerange)
        set(gca,'xlim',timerange);
    end
elseif plottype == 4 %ftimage
    pmax = max(max(max(lfpft(:,:,1:nfreq))));
    [nr,nc] = Nsubplots(length(xvals));
    for k = 1:length(xvals)
        subplot(nr,nc,k);
        imagesc(log(squeeze(lfpft(k,:,1:nfreq))));
        title(sprintf('%.2f',xvals(k)));
        caxis([0 log(pmax)]);
    end
elseif plottype == 5  % stacked time
    subplot(1,1,1);
    hold off;
    for j = 1:size(LFP,2)
        for k = 1:length(xvals)
            plot(lfpt, squeeze(lfpm(k,j,:))+j.*stackoff,'color',colors{j});
            hold on;
        end

    end
elseif plottype == 8  %
    subplot(1,1,1);
    hold off;
        for k = 1:length(xvals)
    a = squeeze(mean(lfpm(k,1:15,:),2));
    b = squeeze(mean(lfpm(k,16:end,:),2));
            plot(lfpt,a ,'r');
            hold on;
            plot(lfpt, b,'b');
            polyfit(a,b,1)
        end

elseif plottype == 6  % stacked FT
    subplot(1,1,1);
    hold off;
    for j = 1:size(LFP,2)
        for k = 1:length(xvals)
            if length(xvals) == 1
                plot(ftfrq(freqs),squeeze(lfpft(k,j,freqs))+j.*0.0,'color',colors{j});
            else
                plot(ftfrq(freqs),squeeze(lfpft(k,j,freqs))+j.*0.0,'color',colors{j});
            end
            hold on;
        end

    end
elseif plottype == 7  %CSD
    subplot(1,1,1);
    csd = diff(mean(LFP,3),2,2);
    imagesc(csd');


elseif plottype ==  9 % plot power vs expt variable
    subplot(1,1,1);
    hold off;
    if calcresp
        fid = find(result.freqs > calcfreq(1) & result.freqs < calcfreq(2)) 
        cresp = sum(lfpft(:,:,:,fid),4);
        imagesc(squeeze(cresp));
        labels = {};
    elseif plotlines
        colors = mycolors;
        hold off;
        nl = 0;
        result.ids = {};

        for ii = 1:prod(env(1:end))
            [cx{:}] = ind2sub(env,ii);
            idx = 1:length(vals{1});
            for k = 1:length(cx)
                if iscellstr(vals{k})
                    id = find(strcmp(uvals{k}{cx{k}},vals{k}));
                else
                    id = find(vals{k} == uvals{k}(cx{k}));
                end
                idx = intersect(idx,id);
            end
            j = cx{1};
            [j,k] = ind2sub([env(1) prod(env(2:end))],ii);
            sc = styles{1+mod(j-1,4)};
            if alln(k,j) > 0
                nl = nl+1;
                result.(extypes{1})(nl) = uvals{1}(cx{1});
                labels{nl} = sprintf('%s=%.2f',et,xvals(j));
                for m = 2:length(cx)
                    result.(extypes{m})(nl) = uvals{m}(cx{m});
                    if iscellstr(uvals{m})
                        labels{nl} = [labels{nl} sprintf(',%s',uvals{m}{cx{m}})];
                    else
                        labels{nl} = [labels{nl} sprintf(',%s=%.2f',extypes{m},uvals{m}(cx{m}))];
                    end
                end
                if isfield(Expt.Trials,'stimid')
                    result.stimid{nl} = unique([Expt.Trials(idx).stimid]);
                end
                labels{nl} = [labels{nl} sprintf(',n=%d',alln(k,j))];
                result.resps(:,:,nl) = squeeze(mean(lfpft(k,j,probes,freqs),3));
                result.n(nl) = alln(k,j);
                result.sid(nl) = j;
                result.ids{nl} = allids{k,j};
                ylines = find(sum(alln') >0);
                yl = find(k ==ylines);
                if length(xvals) == 1
                    plot(result.freqs,squeeze(mean(lfpft(k,j,probes,freqs),3)),'-','color',colors{j});
                else
                    plot(result.freqs,squeeze(mean(lfpft(k,j,probes,freqs),3)),'color',colors{yl},'linestyle',sc);
                end
                hold on;
            end
    end
    legend(labels);
        set(gca,'yscale','log');
        axis('tight');
    else
        resps = squeeze(mean(mean(lfpft,1),2));
        imagesc(resps(:,freqs));
        labels = {};
    end
    lfpm = lfpft;
    result.alln = alln;
    result.labels = labels;
end

result.Header = Expt.Header;
result.Stimvals = Expt.Stimvals;
result.extypes = extypes;