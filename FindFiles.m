function d = FindFiles(monk, varargin)
%FindFiles(monk) finds all usefule data directorys for an animal

dataroot = ['/b/data/' monk];

a = mydir([dataroot '/M*']);
if ~isempty(a)
id = find(CellToMat(regexp({a.filename},'^M[0-8][0-9]+')));
a = a(id);
id = find(~CellToMat(strfind({a.filename},'backup')));
a = a(id);
end
d = a;
a = mydir([dataroot '/G*']);
if ~isempty(a)
id = find(CellToMat(regexp({a.filename},'^G[0-8][0-9]+')));
a = a(id);
id = find(~CellToMat(strfind({a.filename},'backup')));
a = a(id);
d = [d a];
end
d = d(find([d.isdir]));