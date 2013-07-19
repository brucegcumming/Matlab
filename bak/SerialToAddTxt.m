function result = SerialToAddTxt(sfile, ofile, field)


fid = fopen(sfile,'r');
Header.name = sfile;
if fid < 1
    result = [];
    fprintf('Can''t read %s\n',sfile);
    return;
end
tic;
result.filename = sfile;



maxlen = 12000;
a = textscan(fid,'%s','delimiter','\n','bufsize',maxlen);
fclose(fid);
txt = char(a{1});

xid = strmatch(field,txt);
tid = strmatch('O ',txt);

for j = 1:length(xid)
    id = find(tid < xid(j));
    if ~isempty(id)
        t(j) = sscanf(txt(tid(id(end)),:),'O %*d %d');
    end
end
result.t = t;
result.f = txt(xid);
fid = fopen(ofile,'w')
for j = 1:length(t)
    fprintf(fid,'%d %s\n',t(j),deblank(txt(xid(j),:)));
end
fclose(fid);



