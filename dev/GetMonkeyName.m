function [monk, monkeyname, xdir, xname] = GetMonkeyName(name, varargin)
%[monk, monkyename dirsuffix] = GetMonkeyName(name)
%

monk = [];
monkeyname = [];
xname = [];
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
    fprintf('Can''t find monkey name is %s\n',name)
    monk = 'nby';
    monkeyname = 'none';    
end



if regexp(name,['/' monk '/[MG0-9]*/'])
    xdir = regexprep(name,['.*/' monk '/([MG0-9]*)/.*'],'$1');
elseif regexp(name,['/' monk '/[MG0-9]*$'])
    xdir = regexprep(name,['.*/' monk '/([MG0-9]*)'],'$1');
elseif regexp(name,['/' monk '[MG0-9]*/'])
    xdir = regexprep(name,['.*/' monk '([MG0-9]*)/.*'],'$1');
else
    xdir = regexprep(name,['.*/' monk '/(.*)/.*'],'$1');
end
if fullname
    monk = xname;
end
