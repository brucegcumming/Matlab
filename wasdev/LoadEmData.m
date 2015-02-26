
function [cExpt, details]= LoadEmData(cExpt, varargin)
%
%
%  [Expt, details] = LoadEmData(Expt) finds the matching .em matlab file, read in
% the eye position data, and adds these to each trial in Expt
%
% Adds a field EyeData to Each Trials
% Eeydata(:,1) = LH, 2 = RH, 3 = LV 4 = RV
global bgcfileprefix;

  name = strrep(cExpt.Header.Name,'\','/');
  if isfield(cExpt.Header,'Spike2Version')
      [cExpt, details] = LoadSpike2EMData(cExpt, varargin{:});
      return;
  else
  idx = regexp(name,'\.c[0-9]\.');
  emname = strrep(name,name(idx:idx+3),'.em.');
  end
  ok = 1;
  if ~exist(emname,'file')
      emname = [bgcfileprefix emname];
      if ~exist(emname,'file')
          fprintf('No EM file %s\n',emname);
          ok = 0;
      end
  end
  if ok
      load(emname);
      for j = 1:length(Expt.Trials)
          emlen(j) = length(Expt.Trials(j).Eyevals.lh);
      end
      emlen = min(emlen(find(emlen > 0)));
      cExpt.Header.emlen = emlen;
      for j = 1:length(Expt.Trials)
          if Expt.Trials(j).Start == cExpt.Trials(j).Start
              cExpt.Trials(j).EyeData = [];
              cExpt.Trials(j).EyeData(:,1) = Expt.Trials(j).Eyevals.lh(1:emlen);
              cExpt.Trials(j).EyeData(:,2) = Expt.Trials(j).Eyevals.rh(1:emlen);
              cExpt.Trials(j).EyeData(:,3) = Expt.Trials(j).Eyevals.lv(1:emlen);
              cExpt.Trials(j).EyeData(:,4) = Expt.Trials(j).Eyevals.rv(1:emlen);
              cExpt.Trials(j).goodem = 1;
          else
            cExpt.Trials(j).goodem = 0;
          end
      end
  end

  
  
function [cExpt, details] = LoadSpike2EMData(cExpt, varargin)

selectile=5; 
name = strrep(cExpt.Header.Name,'\','/');
emname = strrep(name,'.mat','.em.mat');
if isunix && strncmp(name,'Z:',2)
    emname = strrep(emname,'Z:','/bgc');
elseif isunix && (strncmp(name,'/data',5) || strncmp(name,'/smr',4))
    emname = ['/b' emname];
end
fixid = 0;
reread = 0;
sacargs ={};
suffix = 0;
addsoft = 0;

usetrials = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'addsoft',5)
        addsoft = 1;
    elseif strncmpi(varargin{j},'name',3)
        j = j+1;
        emname = varargin{j};
    elseif strncmpi(varargin{j},'reread',3)
        reread = 1;
    elseif strncmpi(varargin{j},'datdir',6)
        j = j+1;
        datdir = varargin{j};
        emname = strrep(emname,fileparts(emname),datdir);        
    elseif strncmpi(varargin{j},'rmonoc',5)
        sacargs = {sacargs{:} varargin{j}};
    elseif strncmpi(varargin{j},'lmonoc',5)
        sacargs = {sacargs{:} varargin{j}};
    elseif strncmpi(varargin{j},'selectile',9)
        j = j+1;
        selectile = varargin{j};  %set length of record based on this percentile
    elseif strncmpi(varargin{j},'suffix',5)
        j = j+1;
        suffix = varargin{j};
    elseif strncmpi(varargin{j},'trials',5)
        j = j+1;
        usetrials = find(ismember([cExpt.Trials.Trial],varargin{j}));
    elseif strncmpi(varargin{j},'zeropad',5)
        sacargs = {sacargs{:} varargin{j}};
    end
    j = j+1;
end

if isfield(cExpt.Header,'loadname')
    prefix = fileparts(cExpt.Header.loadname);
    [a,b,c] = fileparts(emname);
    newname = [prefix '/' b c];
    if isdir(prefix) && exist(newname)
        emname = newname;
    end
end


if isfield(cExpt.Header,'bysuffix') && cExpt.Header.bysuffix == 1    
    cExpt.Header.bysuffix =2;
    if ~isfield(cExpt.Header,'errs')
        cExpt.Header.errs = {};
    end
    tid = [cExpt.Trials.Trial];
    id = find(diff(tid) < 0);
    if ~isempty(id) 
        if length(id) == length(cExpt.Header.BlockStart) -1; %easy fix
            cExpt.Header.errs{end+1} = sprintf('Fixed Trials.Trial');
            cExpt = FixExpt(cExpt,'Trials');
        else
            cExpt.Header.errs{end+1} = sprintf('Trials.Trial is not consecutive - needs fixing');
        end
        cprintf('error','%s\n',cExpt.Header.errs{end});
    end
    if isfield(cExpt.Header,'BlockStart')
        for j = 1:length(cExpt.Header.BlockStart)
            if j == length(cExpt.Header.BlockStart)
                tid = cExpt.Header.BlockStart(j):cExpt.Trials(end).Trial;
            else
                tid = cExpt.Header.BlockStart(j):cExpt.Header.BlockStart(j+1)-1;
            end
            if isfield(cExpt.Header,'suffixes')
                suffix = cExpt.Header.suffixes(j);
            else
                suffix = cExpt.Header.Combineids(j)
            end
            [cExpt, details{j}] = LoadSpike2EMData(cExpt, varargin{:},...
                'trials',tid,'suffix',suffix);
            for k = 1:length(details{j}.errs)
                cExpt.Header.errs{end+1} = details{j}.errs{k};
            end                
        end
    else
        [cExpt, details] = LoadSpike2EMData(cExpt, varargin{:});
    end
    cExpt.Header.bysuffix =1;
    cExpt = CheckEyeData(cExpt, sacargs{:});
    cExpt = CheckSaccades(cExpt,sacargs{:});

    return;
end
details.errs = {};

if suffix
    emname = regexprep(emname,'\.[0-9]*\.em.mat',['.' num2str(suffix) '.em.mat']);
end

ok = 1;
if ~exist(emname,'file') & exist('bgcfileprefix','var')
    emname = [bgcfileprefix emname];
    if ~exist(emname,'file')
        fprintf('No EM file %s\n',emname);
        ok = 0;
    end
end

if ~exist(emname,'file') && isfield(cExpt.Header,'loadname')
    if isfield(cExpt.Header,'bysuffix') && cExpt.Header.bysuffix == 2
        a = fileparts(cExpt.Header.loadname);
        [b,c] = fileparts(cExpt.Header.loadname);
        c = regexprep(c,'\..*','');
        emname = [a '/' c '.' num2str(suffix) '.em.mat'];
    else
        emname = strrep(cExpt.Header.loadname,'.mat','.em.mat');
    end
end
details.emname = emname;
if isfield(cExpt.Trials,'EyeData') 
    if ~isempty(usetrials)
    elseif reread
        cExpt.Trials = rmfield(cExpt.Trials,'EyeData');
    else
        return;
    end
end


if isempty(usetrials)
    usetrials = 1:length(cExpt.Trials);
end

if ~exist(emname,'file')
    details.errs{end+1} = sprintf('Cannot find %s, needed for Trials %d - %d Id %d - %d\n',...
        emname,usetrials(1),usetrials(end),cExpt.Trials(usetrials(1)).id,cExpt.Trials(usetrials(end)).id);
    [cExpt.Trials(usetrials).goodem] = deal(0);
    [cExpt.Trials(usetrials).emstarti] = deal(NaN);
    cprintf('red','%s\n', details.errs{end})
    return;
end


if ok
    load(emname);

%emtpy start values mess up [Expt.Trials.Start] below
for j = 1:length(Expt.Trials)
    if isempty(Expt.Trials(j).Start)
        Expt.Trials(j).Start = NaN;
    end
    if ~isfield(Expt.Trials,'ftime') || isempty(Expt.Trials(j).ftime)
        Expt.Trials(j).ftime = NaN;
    end
    if isfield(Expt.Trials,'Eyevals')
        V = Expt.Trials(j).Eyevals;
        ls(1) = max(find(~isnan(V.rh)));
        ls(2) = max(find(~isnan(V.lh)));
        ls(3) = max(find(~isnan(V.rv)));
        ls(4) = max(find(~isnan(V.lv)));
        ls = min(ls);
        V.rh = V.rh(1:ls);
        V.lh = V.lh(1:ls);
        V.rv = V.rv(1:ls);
        V.lv = V.lv(1:ls);
        Expt.Trials(j).Eyevals = V;
    end
end


    
    %
% calculation of emlen should also look at start - ftime, so that traces
% are alinged on stimulus on. 
for j = usetrials
    [diffs(j), trial(j)] = min(abs(cExpt.Trials(j).Start(1)-[Expt.Trials.Start]));
    if isfield(Expt.Trials,'Eyevals')
        if isempty(Expt.Trials(trial(j)).Eyevals)
            emlen(j) = 0;
        else
        emlen(j) = min([length(Expt.Trials(trial(j)).Eyevals.rh)...
           length(Expt.Trials(trial(j)).Eyevals.lh) ...
           length(Expt.Trials(trial(j)).Eyevals.rv) ...
           length(Expt.Trials(trial(j)).Eyevals.lv)]);
        end
    else
        emlen(j) = min([length(Expt.Trials(trial(j)).rh)...
            length(Expt.Trials(trial(j)).lh)...
            length(Expt.Trials(trial(j)).rv)...
            length(Expt.Trials(trial(j)).lv)]);
    end
end
if isfield(Expt.Trials,'id') && length(unique([Expt.Trials.id]))
    fixid = 1;
end

%channels can differ in lenght by 1 sample, so take one less that shortest
%for lh
emlens = emlen;
fprintf('%s: Max Start mismatch %.1fms\n',emname,max(diffs./10));
id = find(emlen > 0 & diffs < 500);
emlen = floor(prctile(emlen(id),selectile)) -1;
if isempty(emlen)
    cprintf('red','No Useable Trials in %s\n',emname);
    return;
end
cExpt.Header.emlen = emlen;
%CRsamplerate (ibherited from BW was samples per tick)
cExpt.Header.CRrates = Expt.Header.CRrates;
if Expt.Header.CRrates(1) > 0
    cExpt.Header.CRsamplerate = 1/(10000 * Expt.Header.CRrates(1));
else
    cExpt.Header.CRsamplerate = 597.01./10000; %% rate for uStim files - where eye pos is in wrong channels
end

pre = [Expt.Trials.Start] - [Expt.Trials.ftime];
if isfield(Expt.Header,'preperiod')
prepts = Expt.Header.preperiod./cExpt.Header.CRsamplerate;
else
prepts = 100;
end
cExpt.Header.emtimes = ([1:emlen] .* 1./cExpt.Header.CRsamplerate) - mean(pre(find(~isnan(pre))));
minprepts = floor(min(pre) .* cExpt.Header.CRsamplerate);
if ~isfield(Expt.Trials,'Eyevals')
    for j = 1:length(Expt.Trials)
        Expt.Trials(j).Eyevals.lh = Expt.Trials(j).lh;
        Expt.Trials(j).Eyevals.rh = Expt.Trials(j).rh;
        Expt.Trials(j).Eyevals.rv = Expt.Trials(j).rv;
        Expt.Trials(j).Eyevals.lv = Expt.Trials(j).lv;
    end
end

softoff = [0 0 0 0];
if ~isfield(Expt.Trials,'softoff')
    addsoft = 0;
end
if isfield(Expt.Trials,'Eyevals')
    details.allpretimes = pre;
    pre = [];
    for j = usetrials
        k = trial(j);
        pre(j) = cExpt.Trials(j).Start(1) - Expt.Trials(k).ftime;
    end
    details.pretimes = pre;    
    pretime = min(pre(pre > 0));
    minprepts = floor(pretime .* cExpt.Header.CRsamplerate);
    if minprepts < 1
        cprintf('red','Missing preperiod\n');
    end
    for j = usetrials
        k = trial(j);
        prept(j) = floor((cExpt.Trials(j).Start(1) - Expt.Trials(k).ftime) .* cExpt.Header.CRsamplerate);
    end
    minprepts = floor(prctile(prept(usetrials),5))-1;
    if minprepts < 1
        cprintf('red','Missing preperiod\n');
    end
    for j = usetrials
        starti(j) = prept(j)-minprepts+1;
        if starti(j) <= 0
            lens(j) = emlens(j);
        else
            lens(j) = emlens(j) - starti(j);
        end
    end
    lens(lens<0) = 0;
    maxerr = max(emlen-lens(pre>0));
    emlen = floor(prctile(lens(usetrials),95))-1; 
%emlen is target length total     
    %prept(j)+emlen-minprepts = emlens(j);
%    emlen = emlen-maxerr;
    for j = usetrials
        k = trial(j);
        if addsoft && length(Expt.Trials(k).softoff) >= 4
            softoff = Expt.Trials(k).softoff;
            cExpt.Trials(k).softoff = softoff; 
        end
        if abs(Expt.Trials(k).Start(1) - cExpt.Trials(j).Start(1)) < 500 && starti(j) > -emlen
            cExpt.Trials(j).EyeData = [];
            st = starti(j);
            dst = 1;
            if st <= 0
                dst = 1-st;
                st = 1;
            end
            addpts = emlen-(lens(j)-starti(j));
            if isfield(Expt.Trials(k).Eyevals,'lh')
                if addpts > 0
                    Expt.Trials(k).Eyevals.lh(emlens(j)+1:emlens(j)+addpts) = NaN;
                    Expt.Trials(k).Eyevals.rh(emlens(j)+1:emlens(j)+addpts) = NaN;
                    Expt.Trials(k).Eyevals.lv(emlens(j)+1:emlens(j)+addpts) = NaN;
                    Expt.Trials(k).Eyevals.rv(emlens(j)+1:emlens(j)+addpts) = NaN;
                end
                
                cExpt.Trials(j).EyeData(dst:emlen,1) = Expt.Trials(k).Eyevals.lh(st:emlen+st-dst)+softoff(1);
            cExpt.Trials(j).EyeData(dst:emlen,2) = Expt.Trials(k).Eyevals.rh(st:emlen+st-dst)+softoff(2);
            cExpt.Trials(j).EyeData(dst:emlen,3) = Expt.Trials(k).Eyevals.lv(st:emlen+st-dst)+softoff(3);
            cExpt.Trials(j).EyeData(dst:emlen,4) = Expt.Trials(k).Eyevals.rv(st:emlen+st-dst)+softoff(4);
            if dst > prepts;
                cExpt.Trials(j).EyeData(1:dst-1,:) = NaN;
            end
            else
                fprintf('Missing data in Trial %d\n',k);
            end
            cExpt.Trials(j).goodem = 1;
            cExpt.Trials(j).ftime = Expt.Trials(k).ftime;
            cExpt.Trials(j).emstarti = dst;
            if dst > prepts
                cprintf('red','Missing samples in Trial %d (id%d)\n',cExpt.Trials(j).Trial,cExpt.Trials(j).id);
            end
            if dst > 50
                fprintf('Missing %d samples for Trial %d\n',dst,j);
            end
            Expt.Trials(k).id = cExpt.Trials(j).id;
        else
            cExpt.Trials(j).EyeData(1:emlen,1) = zeros(emlen,1);
            cExpt.Trials(j).EyeData(1:emlen,2) = zeros(emlen,1);
            cExpt.Trials(j).EyeData(1:emlen,3) = zeros(emlen,1);
            cExpt.Trials(j).EyeData(1:emlen,4) = zeros(emlen,1);
            cExpt.Trials(j).goodem = 0;
            cExpt.Trials(j).emstarti = NaN;
        end
        cExpt.Trials(j).emdelay = Expt.Trials(k).Start - cExpt.Trials(j).Start(1);
        szs(j,:) = size(cExpt.Trials(j).EyeData);
    end
else
    for j = usetrials
        k = trial(j);
        if abs(Expt.Trials(k).Start - cExpt.Trials(j).Start) < 300
            cExpt.Trials(j).EyeData = [];
            cExpt.Trials(j).EyeData(1:emlen,1) = Expt.Trials(k).lh(1:emlen)+softoff(1);
            cExpt.Trials(j).EyeData(1:emlen,2) = Expt.Trials(k).rh(1:emlen)+softoff(2);
            cExpt.Trials(j).EyeData(1:emlen,3) = Expt.Trials(k).lv(1:emlen)+softoff(3);
            cExpt.Trials(j).EyeData(1:emlen,4) = Expt.Trials(k).rv(1:emlen)+softoff(4);
            cExpt.Trials(j).ftime = Expt.Trials(k).ftime;
            cExpt.Trials(j).goodem = 1;
        else
            cExpt.Trials(j).EyeData(1:emlen,1) = zeros(emlen,1);
            cExpt.Trials(j).EyeData(1:emlen,2) = zeros(emlen,1);
            cExpt.Trials(j).EyeData(1:emlen,3) = zeros(emlen,1);
            cExpt.Trials(j).EyeData(1:emlen,4) = zeros(emlen,1);
            cExpt.Trials(j).goodem = 0;
        end
    end
end
end
cExpt.Header.emlen = emlen;
cExpt.Header.emtimes = ([1:emlen] .* 1./cExpt.Header.CRsamplerate) - pretime;
cExpt.Header.emChannels = {'lh' 'rh' 'lv' 'rv'}; 
if suffix == 0
evar = mean(var(cat(3,cExpt.Trials(usetrials).EyeData)),3);
if length(sacargs)
    cExpt = CheckSaccades(cExpt,sacargs{:});
elseif (evar(1) > evar(2) * 50 && evar(3) > evar(4) * 10) || (evar(1) > evar(2) * 10 && evar(3) > evar(4) * 50)
    cExpt = CheckSaccades(cExpt,'lmonoc');
elseif evar(2) > evar(1) * 50 && evar(4) > evar(3) * 50
    cExpt = CheckSaccades(cExpt,'rmonoc');
else
    cExpt = CheckSaccades(cExpt);
end
end

if fixid %% Save .em file with Trial.id set correctly
    save(emname,'Expt');
end