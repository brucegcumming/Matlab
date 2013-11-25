function name = GetFullVName(X)

if iscell(X) 
    
elseif ischar(X)
    if strfind(X,'ClusterTimes')
        name = regexprep(X,'ClusterTimes.*','FullV.mat');
    end
end