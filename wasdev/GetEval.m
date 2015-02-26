function [tf, tflist] = GetEval(Expt,type,flag)

% [value, list] =  GetEval(Expt,type,flag)
%
% 'value' returns a stimulus property of type 'type' for an experiment
% This is simple if the property is the same throughout an
% experiment. If more than one value occurs, the mode, median or
% mean value, depending on the stringe flag ('mode', 'median', or 'mean').
% 
% list returns a list of the unique values of the parameter. 
%
% Example, [sf, sflist] = GetEval(Expt,'sf','mean')
% gives the mean value of the spatial frequency, and a list of the
% different vales of SF contained in Expt. 

tnum = 0;
getmode = 0;
if(nargin > 2)
    if isnumeric(flag)
        tnum = flag;
    elseif strmatch(flag,'mode')
    getmode = 1;
  elseif strmatch(flag,'median')
    getmode = 2;
  elseif strmatch(flag,'mean')
    getmode = 0;
  end
else
  getmode = 1;
end

stimtypes = {'none' 'gabor' 'rds' 'grating' 'bar' 'circle' 'rect' 'square' 'probe' ...
    '2grating' 'cylinder' 'corrug' 'sqcorrug' 'twobar' 'rls' 'annulus' 'rdssine' 'nsine' 'rlssine' 'radial' 'image' 'chacker'};

if isempty(Expt)
    tf = NaN;
    return;
end
if strcmp(type,'probesep')
    [tf, tflist] = ReadPenSep(Expt.Header);
    return;
elseif strcmpi(type,'probe')
    if isfield(Expt.Header,'probe')
       tf = Expt.Header.probe;
       tflist = Expt.Header.probe;
    else
        tf = GetProbeFromName(GetEval(Expt,'name'));
    end
    return;
elseif strcmpi(type,'name')
    if isfield(Expt.Header,'expname')
        tflist = Expt.Header.expname;
    else
        tflist = '';
    end
    if isfield(Expt.Header,'loadname')
        tf = Expt.Header.loadname;
    elseif isfield(Expt.Header,'name')
        tf = Expt.Header.name;
    elseif isfield(Expt.Header,'Name')
        tf = Expt.Header.Name;
    end
    return;
end

if(isfield(Expt.Trials,type))
  if length(Expt.Trials(1).(type)) > 1
      x = [];
      for j = 1:length(Expt.Trials)
          if ischar(Expt.Trials(j).(type))
              x{j} = Expt.Trials(j).(type);
          else
          x = [x; Expt.Trials(j).(type)];
          end
      end
  else
      x = cat(1,Expt.Trials.(type));
  end
  if ischar(x)
      if strcmp(type,'Bs')
          for j = 1:size(x,1);
              bs(j) = strmatch(x(j,:),stimtypes,'exact')-1;
          end
          x = bs;
      end
  end

  tflist = sort(unique(x(:)));
  if tnum > 0 && tnum <= length(Expt.Trials)
      tf = Expt.Trials(tnum).(type);
  elseif iscellstr(x)
       [a,b] = Counts(x);     
      tf = b{1};
  elseif(getmode == 1)
    tf = mode(x(:));
  elseif(getmode == 2)
    tf = median(x(:));
  else
    tf = mean(x(:));
  end
elseif isfield(Expt,'Stimvals') & isfield(Expt.Stimvals,type)
  tf = eval(['Expt.Stimvals.' type]);
  if isempty(tf)
      tf = NaN;
  end
  if ischar(tf)
      if strcmp(type,'Bs')
            tf = strmatch(tf,stimtypes,'exact')-1;
      end
      tflist{1} = tf;
  else
  tflist(1) = tf;
  end
elseif strcmp(type,'sz')
    if isfield(Expt.Trials,'wi')
      x = cat(1,Expt.Trials.wi);
      tf = mean(x);
    elseif isisfield(Expt,'Stimvals') & isfield(Expt.Stimvals,'wi')
        tf = Expt.Stimvals.wi;
    end
else
    tf = NaN;
    tflist(1) = tf;
end


function [d, dlist] = ReadPenSep(Header)

d = NaN;             
if isfield(Header,'probesep')
    if isempty(Header.probesep)
        d(2) = NaN;
    else
        d(2) = Header.probesep;
    end
end
if isfield(Header,'Peninfo') && isfield(Header.Peninfo,'trode')
    id = strfind(Header.Peninfo.trode,'Contact');
    if length(id)
        x = id(1);
        id = strfind(Header.Peninfo.trode(x:end),' ');
        if strncmp(Header.Peninfo.trode(id(1)+x:end),'CNT',3)
            x = x+id(2);
        else
            x = x+id(1);
        end
        x = sscanf(Header.Peninfo.trode(x:end),'%d');
        d = x;
    end
elseif isfield(Header,'probesep') & Header.probesep > 0
    d = Header.probesep;
end

if isnan(d)
    [dir, name] = fileparts(strrep(Header.Name,'\','/')); %'/' works on Linux and Windows
    id = strmatch(name(1:7),{'lemM010' 'lemM011'  'lemM019' 'lemM020' 'lemM023' 'lemM027' 'lemM029' 'lemM040' 'lemM043' 'lemM051' 'lemM055'...
        'lemM062' 'lemM067' 'lemM072' 'lemM098' 'lemM106' 'lemM107' 'lemM110'});
    tmpps = [100 100 50 50 50 50 50 150 150 150 150 150 150 150 75 75 75 75];
    if id
        d = tmpps(id(1));
    else
        d= NaN;
    end
end


dlist(1) = d(1);