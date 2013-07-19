function [Expt, success]  = LoadBwLFP(Expt,varargin)

j = 1;
profiling = 0;
plotpwr = 0;
success = 0;
getlfp = 1;
filterargs = {};
lfpname = [];
while j <= length(varargin)
    if strncmpi(varargin{j},'cycleav',5)
        filterargs = {filterargs{:} 'cycleavg'};
    elseif strncmpi(varargin{j},'plotf',5)
        plotpwr = 1;
    elseif ischar(varargin{j})
        lfpname = varargin{j};
    end
    
    j = j+1;
end

if isempty(lfpname)
    name = Expt.Header.Name;
    if ~exist(name)
        name = name2path(name);
    end
    idx = regexp(name,'\.c[0-9]\.');
    lfpname = strrep(name,name(idx:idx+3),'.lfp.');
end
  if ~exist(lfpname,'file')
      lfpname = strrep(lfpname,'bgc','bgc/bgc');
  end
  if ~exist(lfpname,'file')
      fprintf('No file %s\n',lfpname);
      success = 0;
  else
      tic;
    load(lfpname);
    success = 1;
    for j = 1:length(LFP.Trials)
      lfplen(j) = length(LFP.Trials(j).LFP);
    end
    lfpl = prctile(lfplen(find(lfplen > 0)),10);
    lfpmin = min(lfplen(find(lfplen > 0)));
    if lfpmin > lfpl * 0.9
        lfplen = min(lfplen(find(lfplen > 0)));
    else
        id = find(lfplen > lfpl * 0.9);
        lfplen = min(lfplen(id));
        LFP.Trials = LFP.Trials(id);
    end
    Expt.Header.lfplen = lfplen;
    LFP.Header.lfplen = lfplen;
    ftfrq = (0:lfplen-1)/(lfplen/(10000 * LFP.Header.CRsamplerate));
%BW CRsample rate is samples per 0.1ms tic, Spike 2 samplerate is
%secs/sample
    Expt.Header.LFPsamplerate = 1/(LFP.Header.CRsamplerate * 10000);
    k = 1;
    powerspec = zeros(lfplen,1);
    for j = 1:length(LFP.Trials)
        if length(LFP.Trials(j).LFP) > 0
            while max(LFP.Trials(j).LFP) > 5000 %not possible
                [a,b] = max(LFP.Trials(j).LFP);
               LFP.Trials(j).LFP(b) = 0;
            end
%            lfpft(:,j) = fft(LFP.Trials(j).LFP(1:lfplen));
            lfpft = fft(LFP.Trials(j).LFP(1:lfplen));
            powerspec = powerspec + abs(lfpft);                
        end
    end
    if profiling
        toc;
    end
    idx = find(ftfrq > 170 & ftfrq < 190);
    [a,b] = max(powerspec(idx));
    if a > 2 * min(powerspec(idx))
        mainsf = ftfrq(idx(b))/3;
    else
        mainsf = 60;
    end

    Expt.Header.LFPsnr = powerspec(idx(b))./mean(powerspec(2:idx(b)-5));
    if plotpwr
        GetFigure('LFP Power');
        hold off;
        plot(ftfrq(2:400),powerspec(2:400));
        title(sprintf('%s Mains/Signal = %.1f',splitpath(Expt.Header.Name),Expt.Header.LFPsnr));
    end
  if max(powerspec) < 10e-10
      success = 0;
      fprintf('No real LFP data in %s\n',splitpath(Expt.Header.Name));
  end
  end
  
  nrnd = 0;

% first set up dummy variables for things that are encoded in
% Optioncode, make sure RespDir is set for all Trials, and put LFP in if
% required. N.B. LFP.Trials can be larger than Expt.Trials if a cluster was
% only defined for some of the list.

%LFP = FilterLFPTrials(LFP);
k = 1;
for j = 1:length(Expt.Trials);
  if findstr(Expt.Trials(j).OptionCode,'+rp')
    Expt.Trials(j).rndphase = 1;
    nrnd = nrnd +1;
  else
    Expt.Trials(j).rndphase = 0;
  end
  if isfield(Expt.Trials,'RespDir') & isempty(Expt.Trials(j).RespDir);
      Expt.Trials(j).RespDir = 0;
  end
  if getlfp
      while LFP.Trials(k).Start(1) < Expt.Trials(j).Start(1)
          k = k+1;
      end
       if LFP.Trials(k).Start(1) == Expt.Trials(j).Start(1) & ~isempty(LFP.Trials(k).LFP)
  %         Expt.Trials(j).rawlfp = LFP.Trials(k).LFP(1:lfplen);
           [Expt.Trials(j).LFP Expt.Trials(j).FTlfp] = FilterLFP(LFP.Trials(k).LFP(1:lfplen),LFP.Header.CRsamplerate,'freq',mainsf,filterargs);
           Expt.Trials(j).lfptime = LFP.Trials(k).Start(1);
           Expt.Trials(j).lfpo = 0;
           k = k+1;
       end
  end
end
return;
  
