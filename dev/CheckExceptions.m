function CheckExceptions(res, label)
%CheckExceptions(res, label)  runs through a result file
% from AllVPcs, RunAllVpcs, RunAllGridFiles, and prints out 
% any results that were exceptions from try/catch

nerr = 0;

if isstruct(res) && isfield(res,'cls')
    CheckExceptions(res.cls);
    return;
end
if iscell(res)
    nc = prod(size(res));
    for j = 1:nc
        newerr = CheckError(res{j},j);
        if newerr > 0
        nerr = nerr+newerr;
        end
    end
end
if nerr > 0 && nargin > 1  %label given = use a popup
    warndlg(sprintf('%d Errors',nerr),[label  'Exceptions']);
end


function iserr = CheckError(R, id)
    iserr = 0;
    
    if iscell(R)
        for j = 1:length(R)
            CheckError(R{j},id);
        end
        return;
    end
if isfield(R,'errstate')
    s = [];
    line = [];
    if isfield(R,'filename')
        s = [s R.filename];
    end
    if isfield(R,'exptid')
        s = [s ' E' num2str(R.exptid)];
    end
    if isfield(R,'exptid')
        s = [s 'P' num2str(R.probe)];
    end
    if isfield(R,'datadir')
        s = [s R.datadir];
    end
    if id > 0
        s = [s num2str(id)];
    end
    if isempty(s)
        s = 'Unknown Filename';
    end
    [a,b] = fileparts(R.errstate.stack(1).file);
    line = sprintf(' :%s line %d',b,R.errstate.stack(1).line);
    mycprintf('errors','%s:%s%s\n',s,R.errstate.message,line);
    iserr = 1;
end
