function Both = CatStruct(To, From)
% Both = CatStruct(To, From) Add struct elements even if fields don't match


ns = length(To);
f = fields(From);
for j = 1:length(From)
    for k = 1:length(f)
        To(ns+j).(f{k}) = From(j).(f{k});
    end
end
Both = To;

