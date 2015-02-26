function [monk, monkeyname, xdir, xname] = GetMonkeyName(name, varargin)
%[monk, monkyename dirsuffix, xname] = GetMonkeyName(name)
%monk is 3 letter nickname, 
%monkeyname is full name
%xname is the path after the data root
%  [GetFilePath('data') monk xname]  should return correct folder for data 
%                                     from any name

monk = [];
monkeyname = [];
xname = [];
xdir = [];
monkeys = {'lem', 'dae', 'duf', 'jbe', 'bgy', 'ruf', 'ppr', 'ica' 'ic'};
monkeydir = {'lem', 'dae', 'dufus', 'jbe', 'bgy', 'rufus',  'ppr', 'ica' 'icarus'};
fullname = 0;

j = 1;
while j <= length(varargin)
    if strcmpi(varargin{j},'expname')
        fullname = 1;
    end
    j = j+1;
end

if iscell(name)
    for j = 1:length(name)
        [monk{j}, monkeyname{j}, xdir{j}, xname{j}] = GetMonkeyName(name{j}, varargin{:});
    end
    return;
end
if isnumeric(name)
    fprintf('GetMonkeyName Requires a charater string or string array\n',name);
    return;
elseif isstruct(name)
    name = GetName(name);Get
end
name = strrep(name, '\', '/');
for j = 1:length(monkeys)
    if isempty(monk) 
        id = findstr(name, monkeys{j});
        if ~isempty(id)
        monk = monkeys{j};
        monkeyname = monkeydir{j};
%xname is the part of the path that  should not depend on the machine type
        if nargout > 3 || fullname
            if id(1) == 1 %filename not dir
                xname = regexprep(name,'\..*','');
            else
                xname = name(id(1):end);
                xname = regexprep(xname,'/$','');
            end
        end
        end
    end
end

if isempty(monk)
    fprintf('GetMonkeyName:: Can''t find monkey name in %s\n',name)
    monk = 'nby';
    monkeyname = 'none';
end



if regexp(name,['/' monk '/[EMG0-9]+/'])
    xdir = regexprep(name,['.*/' monk '/([EMG0-9]*)/.*'],'$1');
elseif regexp(name,['/' monk '/[EMG0-9]*$'])
    xdir = regexprep(name,['.*/' monk '/([MG0-9]*)'],'$1');
elseif regexp(name,['/' monk '[EMG0-9]+/'])
    xdir = regexprep(name,['.*/' monk '([EMG0-9]+)/.*'],'$1');
elseif regexp(name,['/' monk '/[^/]+/[EMG0-9]+'])
    xdir = regexprep(name,['.*/' monk '/.*/([EMG0-9]+).*'],'$1');
elseif regexp(name,['/' monk '/SE/'])
    xdir = regexprep(name,['.*/' monk '/SE/([SEMG0-9]+).*'],'$1');
elseif isempty(strfind(name,'/'))
    xdir = strrep(name,monk,'');
else
    xdir = regexprep(name,['.*/' monk '/(.*)/.*'],'$1');
    if strcmp(xdir,name) %was no match. Prob just dir name
        xdir = regexprep(name,['.*/' monk '/(.*)'],'$1');
    end
    xdir = regexprep(xdir,['^' monk],'');
end
xdir = regexprep(xdir,'/.*/',''); %lem/XXX/M218 -> M281

if fullname
    monk = xname;
end
