function dtfiles = FindDTFiles(name)

dtfiles = {};
ns = name2path(splitpath(name));
n = 1;
a = regexp(name,'[0-9].c[0-9].');
b = regexp(name,'.mat');
if isempty(a) | isempty(b)
    return;
else
    suffix = name(a(1)+5:b(1)-1);
end

a = strrep(ns,suffix,'rds.OXAC');
if exist(a,'file') dtfiles{n} = a; n = n+1;  end

a = strrep(ns,suffix,'rds.ODX');
if exist(a,'file') dtfiles{n} = a;  n = n+1; end

a = strrep(ns,suffix,'rds.DT');
if exist(a,'file') dtfiles{n} = a;  n = n+1; end

a = strrep(ns,suffix,'rls.OXAC');
if exist(a,'file') dtfiles{n} = a; n = n+1; end

a = strrep(ns,suffix,'rls.ODX');
if exist(a,'file') dtfiles{n} = a;  n = n+1; end

a = strrep(ns,suffix,'rds.AC');
if exist(a,'file') dtfiles{n} = a;  n = n+1; end

