function CompareClusters(Ca,Cb)

fa = fieldnames(Ca);
fb = fieldnames(Cb);
f = intersect(fa,fb);
for j = 1:length(f)
    if isstruct(Ca.(f{j}))
    elseif sum(size(Ca.(f{j}))) ~= sum(size(Cb.(f{j})))
    elseif Ca.(f{j}) ~= Cb.(f{j})
        fprintf('%s\n',f{j});
    end
end