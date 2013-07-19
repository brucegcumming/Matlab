function mnames = GetNames(file, idstr, varargin)

%
% GetNames(file, exptid, suffix, suffix,.....
% returns list of existing files that match the named 'file', but
% have different suffices. 
% for example,
%
% GetNames('/bgc/data/dufus/005/duf005.0.rds.ODX.mat','ODX','SF','SFM')
% returns a cell string array containing
%
% '/bgc/data/dufus/005/duf005.0.rds.SF.mat
% '/bgc/data/dufus/005/duf005.0.rds.SFM.mat
%
% provided these files exist.
%
%
% GetNames('/bgc/data/dufus/005/duf005.0.rds.ODX.mat','ODX','SF','SFM','addgrating')
% returns
% '/bgc/data/dufus/005/duf005.0.rds.SF.mat
% '/bgc/data/dufus/005/duf005.0.rds.SFM.mat
% '/bgc/data/dufus/005/duf005.0.grating.SF.mat
% '/bgc/data/dufus/005/duf005.0.grating.SFM.mat

addgrating = 0;
nvar = nargin - 1;
j = 1;
while j < nvar
  if(strcmpi(varargin{j},'addgrating'))
    addgrating = 1;
  else
    suffixes{j} = varargin{j};
  end
  j = j+1;
end

mnames = {};
fcnt = 1;
grfile = [];

if strfind(file,'rds')
  grfile = strrep(file,'rds','grating');
elseif strfind(file,'corrug')
  grfile = strrep(file,'corrug','grating');
end


if addgrating & length(grfile)
  prefixes = {file, grfile};
  maxk = 2;
else
    prefixes = {file};
    maxk = 1;
end



for k = 1:maxk
  for j = 1:length(suffixes)
  tffile = strrep(prefixes{k},idstr,suffixes{j});
  if(exist(tffile,'file'))
    mnames{fcnt} = tffile;
    fcnt = fcnt + 1;
  end
end
end


