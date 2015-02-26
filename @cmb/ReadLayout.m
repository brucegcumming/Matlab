function DATA = ReadLayout(DATA, file, varargin)   
%think this is obsolete now
DATA = ApplyLayout(DATA, DATA.layoutfile);

return;
layout = [];
if ispc && file(1) == '/'
file(1) = '\';
end
fid = fopen(file,'r');
if fid > 0
a = textscan(fid, '%s','delimiter','\n');
txt = a{1};
for j = 1:length(txt)
sp = findstr(txt{j},' ');
if strncmpi(txt{j},'top',3)
layout.top = sscanf(txt{j}(sp(1)+1:end),'%n');
elseif strncmpi(txt{j},'spkxy',3)
layout.spkxy = sscanf(txt{j}(sp(1)+1:end),'%n');
end
end
fclose(fid);
end


