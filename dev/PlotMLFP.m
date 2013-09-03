function lfpm = PlotMLFP(Expt, varargin)

plottype = 1;
nfreq = 50;
checkdim = 0;
minlen = 0;
chscale = 0;
needft = 0;
stackoff = 0.05;
freqs = 1:nfreq;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'image',4)
        plottype = 2;
    elseif strncmpi(varargin{j},'check',4)
        checkdim = 1;
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            minlen = varargin{j};
        end
    elseif strncmpi(varargin{j},'csd',3)
        plottype = 7;
    elseif strncmpi(varargin{j},'scale',4)
        chscale = 20;
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            chscale = varargin{j};
        end
    elseif strncmpi(varargin{j},'split',4)
        plottype = 8;
    elseif strncmpi(varargin{j},'stack',4)
        plottype = 5;
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1; 
            stackoff = varargin{j};
        end
    elseif strncmpi(varargin{j},'ftpwr',4)
        plottype = 9;
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
    end
    j = j+1;
end

if Expt.Header.rc
    return;
end

et = GetEval(Expt,'et');
eb = GetEval(Expt,'e2');

if ~isfield(Expt.Trials,et)
    [Expt.Trials.(et)] = deal(0);
end

if ~isnan(et)
xiv = [Expt.Trials.(et)];
xvals = unique(xiv);
else
[Expt.Trials.dummy] = deal(1);
xvals = 1;
et = 'dummy';
xiv = [Expt.Trials.dummy];
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
        if ismember(plottype,[3 4 6])
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
    len = size(Expt.Trials(1).LFP,1);
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
else
    PWR = LFP;
end

for j = 1:length(xvals)
    idx = find(xiv == xvals(j));
    for k = 1:size(LFP,2);
        lfpm(j,k,:) = squeeze(mean(LFP(:,k,idx),3));
        lfpft(j,k,:) = squeeze(mean(abs(PWR(:,k,idx)),3));
    end
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
elseif plottype == 2
    [nr,nc] = Nsubplots(length(xvals));
    pmax = max(lfpm(:));
    pmin = min(lfpm(:));
    for k = 1:length(xvals)
        subplot(nr,nc,k);
        imagesc(squeeze(lfpm(k,:,:)));
        title(sprintf('%.2f',xvals(k)));
        caxis([pmin pmax]);
    end
elseif plottype == 4
    pmax = max(max(max(lfpft(:,:,1:nfreq))));
    [nr,nc] = NSubplots(length(xvals));
    for k = 1:length(xvals)
        subplot(nr,nc,k);
        imagesc(squeeze(lfpft(k,:,1:nfreq)));
        title(sprintf('%.2f',xvals(k)));
        caxis([0 pmax]);
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
            plot(ftfrq(freqs),squeeze(lfpft(k,j,freqs))+j.*0.0,'color',colors{j});
            hold on;
        end

    end
elseif plottype == 7  %CSD
    subplot(1,1,1);
    csd = diff(mean(LFP,3),2,2);
    imagesc(csd');


elseif plottype ==  9 % plot power vs expt variae
    subplot(1,1,1);
    hold off;
        for k = 1:length(xvals)
            resps(:,k) = squeeze(mean(lfpft(k,:,freqs),3));
        end
        imagesc(resps);

end