function [result, details] = Recombine(name,varargin)
% Recombine(name,varargin) calls combine to reebuild name
% Recombine(...,'newonly')  only calls ocmbine if some ClusterTimes files
% are newer than name
noquit = 0;
newonly = 0;
args = {};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'newonly',6)
        newonly =1;
    elseif strncmpi(varargin{j},'noquit',6)
        noquit =1;
    end
    j = j+1;
end

if isdir(name)
        d = mydir([name '/*Cells.mat']);
    for j = 1:length(d)
        [res{j}, details{j}]  = Recombine(d(j).name,'noquit',varargin{:});
    end
    combine(details{1}.toplevel,'quit');

elseif ischar(name)
    d = dir(name);
    if isempty(d)
        fprintf('Recombine: Cant Read File %s\n',name);
        return;
    end
    fileage = d.datenum;
    [root, file] = fileparts(name);
    if newonly
        d = dir([root '/Expt*ClusterTimes.mat']);
        id = find([d.datenum] > fileage);
        if isempty(id)
            fprintf('%s is up to date\n',name);
            result.needed = 0;
            result.name = name;
            return;
        end
        for j = 1:length(id)
            fprintf('%s is new\n',d(id(j)).name);
        end
        fprintf('Recombining %s\n',name);
    end
    [a,b,c,d] = GetMonkeyName(name);
    exp = strrep(file,[a c '.'],'');
    if strfind(file,'.Cells')
        mode = 'recombinemucellsquick';
    else
        mode = 'recombinemucells';
    end
    exp = strrep(exp,'.Cells.','');
    exp = strrep(exp,'.Cells','');
    exp = strrep(exp,'Cells.','');
    root = BuildFileName(name,'combine');
    if noquit == 0
        args = {args{:} 'quit'};
    end
    [result, details] = combine(root,'quicksuffix','nospikes','noninteractive','recombinenames',exp,mode,args{:});
    for k = 1:length(result)
       result{k}.ratecheck = PlotRateSequence(result{k});
    end
elseif iscellstr(name)
    for j = 1:length(name)
        result{j} = Recombine(name{j}, varargin{:});
    end
elseif iscell(name) %result of e.g. PlotRateSequence
    for j = 1:length(name)
        if isfield(name{j},'allrates') && isfield(name{j},'missing') %a PlotRateSequence check result
            if ~isempty(name{j}.missing)
                result{j} = Recombine(name{j}.name);
            end
        end
    end
elseif isfield(name,'Expt')&& isfield(name,'Spikes') %an allexpt
    Recombine(name.Expt.Header.loadname,varargin{:});

end