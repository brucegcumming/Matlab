function Trials = ReadSpkFile(name, varargin)

Trials = [];

if ischar(name)
    load(name);
    txt = Ch30;
else
    txt = name;
end


idid = strmatch('id',txt.text);
for j = 1:length(idid)
    ids(j) = sscanf(txt.text(idid(j),3:end),'%f');
end
imvid = strmatch('imve',txt.text);

for j = 1:length(imvid)
    id = find(idid < imvid(j));
    Trials.id(j) = ids(id(end));
    [a,b] = sscanf(txt.text(imvid(j),5:end),'%f,%f %f');
    Trials.imve(j) = a(1);
    Trials.imse(j) = a(2);
end