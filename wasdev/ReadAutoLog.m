function ReadAutoLog(name)

fid = fopen(name);
if fid <= 0
    return;
end

tline = fgets(fid);
nl = 1;
while tline > 0
    id = findstr(tline,'-2011');
    if length(id)
       logt(nl) = datenum(tline(id-6:id+14));
    else
        logt(nl) = 0;
    end
    if strncmp(tline,'Start on',8)
        ops(nl) = 1;
    else
        ops(nl) = 0;
    end
    nl = nl+1;
tline = fgets(fid);
end
fclose(fid);

hold off;
id = find(logt > 0);
plot(id,(logt(id)-logt(id(1))).*24);
totaltime = (logt(id(end)) - logt(id(1)))*24;
hold on;
id = find(logt > 0 & ops == 1);
plot(id,(logt(id)-logt(id(1))).*24,'ro');

id = find(ops(2:end) ==1);
d = diff(logt);
buildtime = sum(d(id)) .*24;
fprintf('Build times total %.2f/%.2f',buildtime,totaltime);

