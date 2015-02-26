function F = split(S,c)
%cellstr = split(s,delimiter)  splits a string into a cell array

if nargin == 1
    c = '\s+';
    nc = 1;
else
    nc = length(c);
    if strcmp(c,'\n')
        nc = 1;
    end
end
id = regexp(S,c);
%consider a different loop for nc > 1
%can be a lot slower 
if isempty(id)
    F{1} = S;
elseif nc > 1
    last = 1;
    for j = 1:length(id)
        F{j} = S(last:id(j)-1);
%        F{j} = regexprep(F{j},['^' c],''); %makes it very slow
        last = id(j)+nc;
    end
    if last <= length(S)+nc-1
        j = j+1;
        if last >length(S)
            F{j} = '';
        else
            F{j} = S(last:end);
        end
%        F{j} = regexprep(F{j},['^' c],'');
    end
else
    last = 1;
    for j = 1:length(id)
        F{j} = S(last:id(j)-1);
%        F{j} = regexprep(F{j},['^' c],''); %makes it very slow
        last = id(j)+1;
    end
    if last <= length(S)
        j = j+1;
        if last >length(S)
            F{j} = '';
        else
            F{j} = S(last:end);
        end
%        F{j} = regexprep(F{j},['^' c],'');
    end
end