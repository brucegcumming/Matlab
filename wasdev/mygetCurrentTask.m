function t = mygetCurrentTask(varargin)
%like current task, but returns valid structure when no pool open.
% mygetCurrentTask('print') printts out worker ID if its nonzero
% mygetCurrentTask('num') just returns worker number
print = 0;
numeric = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'numeric',3)
        numeric = 1;
    elseif strncmpi(varargin{j},'print',5)
        print = 1;
    end
    j = j+1;
end

t = getCurrentTask();
if isempty(t)
    t.ID = 0;
elseif print
    fprintf('Worker %d\n',t.ID);
end


if numeric
    t = t.ID;
end