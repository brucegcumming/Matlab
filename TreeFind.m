function [res, sizes, dates, extra] = TreeFind(path, varargin)
% [names, sizes, dates, pathnames] = TreeFind(path, varargin)
% returns a list of files under the tree starting in Path.
%
% TreeFind(path,'name','xxx')   finds files where strfind('xxx',pathname)
% is true
%
%
% TreeFind(path,'function','myfunc')   names a function to be executed. 
% Finds filew with return value >1 finds files where strfind('xxx',pathname)
% TreeFind(path,'function','myfunc',{args})   names a function to be executed. 
% calls myfunc with myfunc(file, args{:});
%
%   ....,'name','.*[0-9][0-9].mat') finds raw spk2 .mat files
name = [];
funcfcn = [];
funcargs = [];
printfiles = 0;
findzero = 1;
found = 0;
newer = 0;
res = {};
sizes = [];
dates = [];
extra = [];
recursive = 1;
verbose = 0;
exclude = {};

j = 1;
while j < nargin
    if strncmpi(varargin{j},'name',3)
        j = j+1;
        name = varargin{j};
    elseif strncmpi(varargin{j},'newer',3)
        j = j+1;
        newer = varargin{j};
    elseif strncmpi(varargin{j},'exclude',3)
        j = j+1;
        exclude{end+1} = varargin{j};
    elseif strncmpi(varargin{j},'nonrecursive',5)
        recursive = 0;
    elseif strncmpi(varargin{j},'function',3)
        j = j+1;
        funcfcn = varargin{j};
        if length(varargin) > j & iscell(varargin{j+1})
            j = j+1;
            funcargs = varargin{j};
        end
    elseif strncmpi(varargin{j},'print',3)
        printfiles = 1;
    elseif strncmpi(varargin{j},'verbose',5)
        verbose = 1;
    end
    j = j+1;
end

d = dir(path);
for j = 1:length(d)
   if d(j).isdir & isempty(strmatch('.',d(j).name)) && isempty(strmatch('..',d(j).name)) && recursive > 0
       if sum(cellstrfind(d(j).name,exclude)) == 0
        if isdir(path)
            subname = [path '/' d(j).name];
        else  %wildcard in path
            subname = [fileparts(path) '/' d(j).name];
        end
        if verbose
            fprintf('Searching %s\n',subname);
        end
        [files, ns, nd, nu] = TreeFind(subname,varargin{:});
        res = {res{:} files{:}};
        sizes = [sizes ns];
        dates = [dates nd];
        extra = [extra nu];
        found = length(res);
       end
   else
       good = 1;
        if ~isempty(name) & isempty(regexp(d(j).name,name))
            good = 0;
        elseif isempty(d(j).datenum)
            good = 0;
        else
            good = 1;
        end
        if newer > 0 && (now-d(j).datenum) > newer
            good = 0;
        end
        pathname = [path '/' d(j).name];
        if good & ~isempty(funcfcn)
            if isempty(funcargs)
            good = eval([funcfcn '(''' pathname ''')']);
            else
                str = [];
                for k = 1:length(funcargs)
                    str = [str ',''' funcargs{k} ''''];
                end
                good = eval([funcfcn '(''' pathname '''' str ')']);
            end
            if isstruct(good)
            elseif isnan(good)
                nu = length(extra)+1;
                extra(nu).unknown = pathname;
                good = 0;
            else
                if good == 0 & findzero
                    good = 1;
                end
            end
        end
        if good
            found = found + 1;
            res{found} = pathname;
            sizes(found) = d(j).bytes;
            dates(found) = d(j).datenum;
            extra(found).isdir = d(j).isdir;
            if printfiles
                fprintf('%s\n',pathname);
            end
        end
    end
end