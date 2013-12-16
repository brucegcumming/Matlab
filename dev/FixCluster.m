function [C, errs] = FixCluster(C, varargin)
%FixCluster(C, varargin) modifids old clustr structs to make sure
%FixCluster(C, DATA) checks fields in DATA that should be in Cluster, like
%chspk
errs = {};
nerr = 0;
DATA = [];

j = 1;
while j <=length(varargin)
    if isfield(varargin{j},'chspk')
        DATA = varargin{j};
    end
    j = j+1;
end

if iscell(C)
    for j = 1:length(C)
        if ~isfield(C{j},'probe')
            C{j}.probe = j;
        end
        [C{j}, err] = FixCluster(C{j},varargin{:});
        errs = {errs{:} err};
    end
    return;
end


if isfield(DATA,'chspk') && ~isfield(C,'chspk')
    C.chspk = DATA.chspk;
end

if ~isfield(C,'cluster')
    nerr = nerr+1;
    errs{nerr} = 'missingcluster';
    C.cluster = 1;
end
if ~isfield(C,'manual')
    nerr = nerr+1;
    errs{nerr} = 'missingmanual';
    C.manual = 0;
end
if ~isfield(C,'auto')
    nerr = nerr+1;
    errs{nerr} = 'missingauto';
    C.auto = 0;
end
if ~isfield(C,'savetime')
    C.savetime = [0 0];
end
if ~isfield(C,'exptno')
    nerr = nerr+1;
    errs{nerr} = 'missingexptno';
    C.exptno = 0;
    if isfield(C,'spkfile')
        id = regexp(C.spkfile,'*p[0-9]*t[0-9]*');
        if ~isempty(id)
        [a,C.exptno] = sscanf(C.spkfile(id(1):end),'p%dt%d');
        end
    end
end
if ~isfield(C,'next')
    C.next = {};
    nerr = nerr+1;
    errs{nerr} = 'missingnext';
elseif ~iscell(C.next)
    nerr = nerr+1;
    errs{nerr} = 'next struct';
    next = C.next;
    C = rmfield(C,'next');
    nx = 1;
    C.next{nx} = rmfields(next,'next');
    while isfield(next,'next')
        nx = nx+1;
        next = next.next;
        C.next{nx} = rmfields(next,'next');
    end
end

