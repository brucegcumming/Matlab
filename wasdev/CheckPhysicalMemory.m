function mem = CheckPhysicalMemory()
%find physical memory in Mb, using 'feature memstats'

if ispc
a = evalc('feature memstats');
id = strfind(a,'Total');
mem = sscanf(a(id(1)+10:end),'%d');
else
    mem = 0;
end