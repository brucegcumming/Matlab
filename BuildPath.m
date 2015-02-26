function path = BuildPath(name, varargin)
%path = BuildPath(name, varargin)
    drives = {'Y:' 'Z:' 'X:' 'C:'};
pathonly = 0;
j = 1;
while j <= length(varargin)
    if strcmpi(varargin{j},'local')
        drives = {'C:' 'Y:' 'Z:' 'X:'};
    elseif strncmpi(varargin{j},'dir',3)
        pathonly = 1;
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
rest = regexprep(rest,['^' mnk '/' num '/'],'');
if pathonly
    rest = [];
end
if ispc 
    needdrive = 0; %if /b/x works on PC, don't need drive. Or if name already has drive
    if needdrive
        for j = 1:length(drives)
            if isdir([drives{j} '/b/data'])
                path = [drives{j} '/b/data/' mnk '/' num '/' rest];
                break;
            end
        end
    else
        path = ['/b/data/' mnk '/' num '/' rest];
    end
else
    path = ['/b/data/' mnk '/' num '/' rest];
end