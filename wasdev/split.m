function F = split(S,c)
%cellstr = split(s,delimiter)  splits a string into a cell array

if nargin == 1
    c = '\s+';
end
id = regexp(S,c);
if isempty(id)
    F{1} = S;
else
    last = 1;
    for j = 1:length(id)
        F{j} = S(last:id(j)-1);
        F{j} = regexprep(F{j},['^' c],'');
        last = id(j)+1;
    end
    if last <= length(S)
        j = j+1;
        F{j} = S(last:end);
        F{j} = regexprep(F{j},['^' c],'');
    end
end