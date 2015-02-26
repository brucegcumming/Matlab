function path = BuildPath(name, varargin)
    drives = {'Y:' 'Z:' 'X:' 'C:'};

j = 1;
while j <= length(varargin)
    if strcmpi(varargin{j},'local')
        drives = {'C:' 'Y:' 'Z:' 'X:'};
    end
    j = j+1;
end
    [mnk, mnkname, num, rest] = GetMonkeyName(name);

d = [mnk '/' num '/'];
if strncmp(num,name,length(num)) % not a path, just a name
    rest = num;
    num = strrep(name,mnk,'');
    num = regexprep(num,'\..*','');
end
if ispc
    for j = 1:length(drives)
        if isdir([drives{j} '/b/data'])
            path = [drives{j} '/b/data/' mnk '/' num '/' rest];
            break;
        end
    end
else
    path = ['/b/data/' mnk '/' num '/' rest];
end