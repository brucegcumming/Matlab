function took = fseektest(name, nloops, dsize)
% fseektest(name, nloop)
%
% opens a read file handle on name, and
% then performs nloop fseek/fread operation
% from random locations

fid = fopen(name,'r');
if fid < 0
    fprintf('Cant open %s\n',name);
    return;
end
ndata = ftell(fid);
ns = ceil(rand(1,nloops) .*ndata-dsize);
ts = now;
for j = 1:nloops
    if j > 1
        fseek(fid, ns(j), 'bof');
    end
    x = fread(fid, dsize,'int16');
end
fclose(fid);
took = mytoc(ts);