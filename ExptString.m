function str = ExptString(E, s)
%str = ExptString(E, s) Return a string descibing property s of Expt

str = [];

if iscell(E)
    for j = 1:length(E)
        str{j} = ExptString(E{j},s);
    end
    return;
end

if iscellstr(s)
    for j = 1:length(s)
        str = [str ExptString(E,s{j})];
    end
    return;
end

x = GetEval(E,s);
str = [s num2str(x)];
