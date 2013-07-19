function CompareExpts(Ea,Eb, varargin)


f = fields(Ea.Stimvals);
for j = 1:length(f);
    if Ea.Stimvals.(f{j}) ~= Eb.Stimvals.(f{j})
        fprintf('%s different\n',f{j});
    end
end
