function Vs = LoadAllFullV(name, varargin)
verbose = 0;
Vs = {};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'verbose',5)
        verbose = 1;
    end
    j = j+1;
end

ts = now;
d = mydir([name '/*FullV.mat']);
bytes = 0;
for j = 1:length(d)
    if verbose
        fprintf('Loading %s\n',d(j).name);
    end
    Vs{j} = LoadFullV(d(j).name,'noconvert');
    bytes = bytes + d(j).bytes;
end
t = mytoc(ts);
fprintf('%d files took %.1f sec = %.2f Mb/sec \n',length(d),t,bytes./(t.*2^20));