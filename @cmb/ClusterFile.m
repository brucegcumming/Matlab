function cfile = ClusterFile(DATA,varargin)
getonline = 0;
getauto = 0;
allprobes = 0;
spkcodefile = 0;
probe = DATA.probe;
j = 1;
while j <= length(varargin)
if strncmpi(varargin{j},'auto',4) %Read Clusters from Online mat file
getauto = 1;
elseif strncmpi(varargin{j},'allprobes',4) %Read Clusters from Online mat file
allprobes = 1;
elseif strncmpi(varargin{j},'getonline',4) %Read Clusters from Online mat file
getonline = 1;
elseif strncmpi(varargin{j},'codes',4) %Save codes/times
spkcodefile = 1;
elseif strncmpi(varargin{j},'probe',4) %Read Clusters from Online mat file
j = j+1;
probe = varargin{j};
end
j = j+1;
end
if DATA.state.online == 0 %%not online data
if getonline
[a,b] = splitpath(DATA.datafilename);
if length(b)
cfile = [b '/OnlineClusters.mat'];
else
cfile = 'OnlineClusters.mat';
end
elseif getauto
cfile = strrep(DATA.datafilename,'.mat','.autocl.mat');
elseif length(DATA.probelist) > 1 && allprobes
cfile = strrep(DATA.datafilename,'.mat','.allcl.mat');
elseif spkcodefile 
cfile = strrep(DATA.datafilename,'.mat','codes.mat');
elseif length(DATA.probelist) > 1 
cfile = strrep(DATA.datafilename,'.mat',['.p' num2str(probe) 'cl.mat']);
else
cfile = strrep(DATA.datafilename,'.mat','.cl.mat');
end
else
if spkcodefile 
cfile = sprintf('%s/SpikeCodes%d.mat',DATA.datafilename,DATA.Expts{DATA.currentexpt(1)}.Trials(end).id);
else
cfile = [DATA.datafilename '/Clusters.mat'];
end
end

