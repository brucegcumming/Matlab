function FixEdepth(name,varargin)

savefix = 0;
if isstruct(name) && isfield(name,'text')
    Ch30 = name;
else
load(name);
end
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'save',4)
        savefix = 1;
    end
end
Ch30.text(:,8) = ' ';
tid = strmatch('ed',Ch30.text);
eds = sscanf(Ch30.text(tid,3:8)','%f');
de = diff(eds);
id = find(de < -0.15);
de(id) = 0;
de = [0; de];
fixed = cumsum(de)+eds(1);
plot(eds,'o-');
hold on;
plot(fixed,'r');
    for j = 1:length(tid)
Ch30.text(tid(j),1:8) = sprintf('ed%.4f',fixed(j));
    end
if savefix
    clear eds;
    clear de;
    clear id;
    clear tid;
    clear fixed;
    clear savefix;
    save(name);
end