function str = BuildFileName(name, type, varargin) 
%str = BuildFileName(name, type, ...) returns string naming file of type
% types are 'fullv','combine','rffits' 'datadir' 'celllist'
%
str = [];
probe  =0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'probe',5)
        j = j+1;
        probe = varargin{j};
    end
    j = j+1;
end
if iscell(name)
    for j = 1:length(name)
        str{j} = BuildFileName(name{j},type, varargin{:});
    end
    return;
end

eid = GetExptNumber(name);
if isfield(name,'Header')
    name = GetEval(name,'name');
elseif ~ischar(name)
    name = GetName(name);
end

[a,b,c,d] = GetMonkeyName(name);
if strcmp(type,'fullv')
    if probe
        str = [fileparts(name) '/Expt' num2str(eid) '.p' num2str(probe) 'FullV.mat'];
    else
        str = [fileparts(name) '/Expt' num2str(eid) 'FullV.mat'];
    end
elseif strcmp(type,'combine') % make filename to give to combine
    [root, file] = fileparts(name);
    if c(1) == 'M' %laminar probe - use directory
        str = root;
    else
        str = [root '/' a c '.mat'];
    end
elseif strcmp(type,'datadir')
    str = regexprep(name,['/' c '/.*'],['/' c]);
    if isempty(strfind(name,'/'))
        str = ['/b/data/' a '/' c '/'];
    end
elseif strcmp(type,'celllist')
    str = regexprep(name,['/' c '/.*'],['/' c]);
    if isempty(strfind(name,'/'))
        str = ['/b/data/' a '/' c '/'];
    end
    str = [str '/CellList.mat'];
elseif strcmp(type,'rffits')
    if regexp(name,[a '/SE/'])
        dataroot = regexprep(name,['data/' a '/SE/.*' ],['data/' a '/SE/']);
    else
        dataroot = regexprep(c,['data/' a '/.*' ],['data/' a]);
    end
    str = [dataroot '/rffits.mat'];
elseif strcmp(type,'bnc')
    str = regexprep(name,['/' c '/.*'],['/' c]);
elseif strcmp(type,'rffile')
    if regexp(name,[a '/SE/'])
        dataroot = regexprep(name,['data/' a '/SE/.*' ],['data/' a '/SE/']);
    else
        dataroot = regexprep(name,['data/' a '/.*' ],['data/' a]);
    end
    str = [dataroot '/' c '/' a c '.rf.mat'];
end
