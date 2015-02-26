function PrintLogData(DATA, type)
id = find(DATA.log.types == type);
if isempty(id)
return;
end
allexp = unique(cat(1,DATA.log.savexps{id}));
[done, b] = max(DATA.log.savedates(id));
last = id(b);
label = DATA.log.expnames{last};
fprintf('%s %s: Expts ',label,datestr(done));
fprintf('%d ',DATA.log.savexps{last});
x = setdiff(allexp, DATA.log.savexps{last});
if length(x)
fprintf('Other '); 
fprintf('%d ',x);
end
fprintf('\n');



