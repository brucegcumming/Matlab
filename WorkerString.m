function s = WorkerString(varargin)
%Return string with worker number 
j = 1;
while j <= length(varargin)
    j = j+1;
end

t = getCurrentTask();
if isempty(t)
    s = '0';
elseif 
    s = sprintf('%d',t.ID);
end
