function p = Name2Probe(name)
p = 0;
id = regexp(name,'\.p[0-9]*');
if ~isempty(id)
    p = sscanf(name(id(1)+2:end),'%d');
end
