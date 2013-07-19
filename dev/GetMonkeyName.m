function [monk, monkeyname, xdir] = GetMonkeyName(name)
%[monk, monkyename] = GetMonkeyName(name)
%

monk = [];
monkeyname = [];
monkeys = {'lem', 'dae', 'duf', 'jbe', 'bgy', 'ruf', 'ppr', 'ica' 'ic'};
monkeydir = {'lem', 'dae', 'dufus', 'jbe', 'bgy', 'rufus',  'ppr', 'ica' 'icarus'};
    

for j = 1:length(monkeys)
    if isempty(monk) && ~isempty(findstr(name, monkeys{j}))
        monk = monkeys{j};
        monkeyname = monkeydir{j};
    end
end

if isempty(monk)
    fprintf('Can''t find monkey name is %s\n',name)
    monk = 'nby';
    monkeyname = 'none';    
end

if regexp(name,['/' monk '/[MG0-9]*/'])
    xdir = regexprep(name,['.*/' monk '/([MG0-9]*)/.*'],'$1');
elseif regexp(name,['/' monk '[MG0-9]*/'])
    xdir = regexprep(name,['.*/' monk '([MG0-9]*)/.*'],'$1');
else
    xdir = regexprep(name,['.*/' monk '/(.*)/.*'],'$1');
end
