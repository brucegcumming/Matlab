
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

name = strrep(cExpt.Header.Name,'\','/');
emname = strrep(name,'.mat','.em.mat');
if isunix && strncmp(name,'Z:',2)
    emname = strrep(emname,'Z:','/bgc');
elseif isunix && (strncmp(name,'/data',5) || strncmp(name,'/smr',4))
    emname = ['/bgc' emname];
end
fixid = 0;
reread = 0;
sacargs ={};
suffix = 0;

usetrials = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'name',3)
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
    elseif strncmpi(varargin{j},'suffix',5)
        j = j+1;
        suffix = varargin{j};
    elseif strncmpi(varargin{j},'trials',5)
        j = j+1;
        usetrials = find(ismember([cExpt.Trials.Trial],varargin{j}));
    end
    j = j+1;
end


if isfield(cExpt.Header,'bysuffix') && cExpt.Header.bysuffix == 1    
    cExpt.Header.bysuffix =2;
    for j = 1:length(cExpt.Header.BlockStart)
        if j == length(cExpt.Header.BlockStart)
            tid = cExpt.Header.BlockStart(j):cExpt.Trials(end).Trial;
        else
            tid = cExpt.Header.BlockStart(j):cExpt.Header.BlockStart(j+1)-1;
        end
        [cExpt, details{j}] = LoadSpike2EMData(cExpt, varargin{:},...
            'trials',tid,'suffix',cExpt.Header.Combineids(j));
    end
    cExpt.Header.bysuffix =1;
    cExpt = CheckEyeData(cExpt);
    cExpt = CheckSaccades(cExpt,sacargs{:});

    return;
end

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
if isfield(cExpt.Trials,'EyeData') 
    if ~isempty(usetrials)
    elseif reread
        cExpt.Trials = rmfield(cExpt.Trials,'EyeData');
    else
        return;
    end
end

if ~exist(emname,'file')
    return;
end

if isempty(usetrials)
    usetrials = 1:length(cExpt.Trials);
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
emlen = floor(prctile(emlen(id),95)) -1;
if isempty(emlen)
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
cExpt.Header.emtimes = ([1:emlen] .* 1./cExpt.Header.CRsamplerate) - mean(pre(find(~isnan(pre))));
minprepts = floor(min(pre) .* cExpt.Header.CRsamplerate);

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
    for j = usetrials
        k = trial(j);
        prept(j) = floor((cExpt.Trials(j).Start(1) - Expt.Trials(k).ftime) .* cExpt.Header.CRsamplerate);
    end
    minprepts = floor(prctile(prept,5))-1;
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
    emlen = floor(prctile(lens(usetrials),5))-1; 
%emlen is target length total     
    %prept(j)+emlen-minprepts = emlens(j);
%    emlen = emlen-maxerr;
    for j = usetrials
        k = trial(j);
        if abs(Expt.Trials(k).Start(1) - cExpt.Trials(j).Start(1)) < 500
            cExpt.Trials(j).EyeData = [];
            st = starti(j);
            dst = 1;
            if st <= 0
                dst = 1-st;
                st = 1;
            end
            addpts = emlen-(lens(j)-starti(j)+1);
            if isfield(Expt.Trials(k).Eyevals,'lh')
                if addpts > 0
                    Expt.Trials(k).Eyevals.lh(emlens(j):emlens(j)+addpts) = NaN;
                    Expt.Trials(k).Eyevals.rh(emlens(j):emlens(j)+addpts) = NaN;
                    Expt.Trials(k).Eyevals.lv(emlens(j):emlens(j)+addpts) = NaN;
                    Expt.Trials(k).Eyevals.rv(emlens(j):emlens(j)+addpts) = NaN;
                end
                
                cExpt.Trials(j).EyeData(dst:emlen,1) = Expt.Trials(k).Eyevals.lh(st:emlen+st-dst);
            cExpt.Trials(j).EyeData(dst:emlen,2) = Expt.Trials(k).Eyevals.rh(st:emlen+st-dst);
            cExpt.Trials(j).EyeData(dst:emlen,3) = Expt.Trials(k).Eyevals.lv(st:emlen+st-dst);
            cExpt.Trials(j).EyeData(dst:emlen,4) = Expt.Trials(k).Eyevals.rv(st:emlen+st-dst);
            else
                fprintf('Missing data in Trial %d\n',k);
            end
            cExpt.Trials(j).goodem = 1;
            cExpt.Trials(j).emstarti = minprepts;
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
    end
else
    for j = usetrials
        k = trial(j);
        if abs(Expt.Trials(k).Start - cExpt.Trials(j).Start) < 300
            cExpt.Trials(j).EyeData = [];
            cExpt.Trials(j).EyeData(:,1) = Expt.Trials(k).lh(1:emlen);
            cExpt.Trials(j).EyeData(:,2) = Expt.Trials(k).rh(1:emlen);
            cExpt.Trials(j).EyeData(:,3) = Expt.Trials(k).lv(1:emlen);
            cExpt.Trials(j).EyeData(:,4) = Expt.Trials(k).rv(1:emlen);
            cExpt.Trials(j).goodem = 1;
        else
            cExpt.Trials(j).EyeData(:,1) = zeros(emlen,1);
            cExpt.Trials(j).EyeData(:,2) = zeros(emlen,1);
            cExpt.Trials(j).EyeData(:,3) = zeros(emlen,1);
            cExpt.Trials(j).EyeData(:,4) = zeros(emlen,1);
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