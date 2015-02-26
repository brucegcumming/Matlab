function [name, details] = GetName(C, varargin)
% Getname(X) returns the cell name associated with string/structre X
% e.g. GetName(Expt) will return lemM001 if thats the correct Expt session
%Getname(X,'path') returns a full pathname
%Getname(X,'folder') returns a the folder the file is in
%if X is a cell array, builds a cell array of names.  If
% this has one unique value, returns this value. Otherwise returns the cell
% array

fullpath = 0;
nametype = '';
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'path',4)
        fullpath = 1;
    elseif strncmpi(varargin{j},'folder',6)
        fullpath = 2;
    elseif strncmpi(varargin{j},'lfp',3)
        nametype = 'lfpfile';
    end
    j = j+1;
end

details = [];
name = '';
if isstruct(C)
    if isfield(C,'Expt') %Expt
        [name, details] = GetName(C.Expt, varargin{:});
    elseif isfield(C,'Header') %Expt
        ename = GetEval(C,'name');
        if fullpath
            name = ename;
        else
            [a,b,c,d] = GetMonkeyName(ename);
            name = [a c];
        end
    elseif isfield(C,'dirname') %RF fit
        [a,b,c,d] = GetMonkeyName(C.dirname);
        name = [a c];
    elseif isfield(C,'V') %FullV file
        if isfield(C,'loadname')
            [a,b,c,d] = GetMonkeyName(C.loadname);
        elseif isfield(C,'name')
            [a,b,c,d] = GetMonkeyName(C.name);
        end
        name = [a c];
        if fullpath
            name = [GetFilePath('data') '/' d];
        end
    elseif isfield(C,'loadname') %could be many thing. including RF fix list
        [a,b,c,d] = GetMonkeyName(C.loadname);
        name = [a c];
    elseif isfield(C,'name') %could be many thing. including RF fix list
        [a,b,c,d] = GetMonkeyName(C.name);
        name = [a c];
    end
elseif ischar(C)
        [a,b,c,d] = GetMonkeyName(C);
        name = [a c];
elseif iscell(C)
    for j = 1:length(C)
        names{j} = GetName(C{j},varargin{:});
    end
    uname = unique(names);
    if length(uname) ==1
        name = uname{1};
        details.names = names;
    else
        name = names;
    end
    return;
end
if fullpath ==2
    name = fileparts(name);
end


if strcmp(nametype,'lfpfile')
    eid = GetExptNumber(C);
    name = [name '.' num2str(eid) '.lfp.mat'];
end


