function cplists(caller,b)

lfig = get(caller,'parent');
dat = get(lfig,'UserData');
it = findobj(lfig,'Tag','SuffList');
for j = 1:length(it)
go = get(it(j),'value');
if go
fprintf('Adding %s to %s\n',dat.files{j},dat.lists{j});
AddToList(dat.lists{j},dat.files{j});
end
end
close(lfig);

