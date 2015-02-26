function Trials = ReadDxseq(name)
%get dx: lines from serial output and match up with id numbers

if ischar(name)
    lines = scanlines(name);
elseif iscellstr(name)
    lines = name;
end

rfid = find(strncmp('RF dx', lines,5));
dxid = find(strncmp('dx:', lines,3));
ceid = find(strncmp('ce:', lines,3));
idid = find(strncmp('id', lines,2));
exid = find(strncmp('exvals', lines,6) & ~strncmp('exvals ', lines,7));

lastid = 0;

for j = 1:length(rfid)
    id = find(dxid < rfid(j) & dxid > lastid);
    if ~isempty(id)
        dxline = dxid(id(end));
        nc = 4;
        a = sscanf(lines{dxline}(nc:end),'%f ');
        Trials(j).dx = a;
        lastid = dxline;
        id = find(dxid < rfid(j));
        line = lines{ceid(id(end))};
        a = sscanf(line(nc:end),'%f ');
        Trials(j).ce = a;
        Trials(j).result = 0;
        
        idline = find(idid < rfid(j));
        a = sscanf(lines{idid(idline(end))},'id%d');
        Trials(j).id = a;

        x = find(exid < rfid(j));
        a = sscanf(lines{exid(x(end))}(7:end),'%f');
        Trials(j).exvals = a;
    else        
        Trials(j).result = -10;
    end
end
id = find([Trials.result] ==0);    
Trials = Trials(id);