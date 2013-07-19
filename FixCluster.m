function [C, errs] = FixCluster(C, varargin)
%FixCluster(C, varargin) modifids old clustr structs to make sure
errs = {};
nerr = 0;


if iscell(C)
    for j = 1:length(C)
        [C{j}, err] = FixCluster(C{j},varargin{:});
        errs = {errs{:} err};
    end
    return;
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

