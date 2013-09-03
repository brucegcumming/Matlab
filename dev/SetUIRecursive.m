function SetUIRecursive(F, S, varargin)


f = fields(S);
FF = get(F);
for j = 1:length(f)
    if isfield(FF,f{j})
        set(F,f{j},S.(f{j}));
    end
end
c = get(F,'children');
for j = 1:length(c)
    SetUIRecursive(c(j),S);
end


