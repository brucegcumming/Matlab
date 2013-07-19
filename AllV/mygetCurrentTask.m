function t = mygetCurrentTask()
%liek current task, but returns valid structure when no pool open.

t = getCurrentTask();
if isempty(t)
    t.ID = 0;
end