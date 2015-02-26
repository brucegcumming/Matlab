function ClearAnnotations(a)
c = get(a,'Children')
for j = 1:length(c)
type = get(c(j),'Type');
if strcmp(type,'line')
delete(c(j));
end
end

