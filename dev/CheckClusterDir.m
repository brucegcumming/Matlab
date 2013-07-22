function DATA = CheckClusterDir(name, varargin)

DATA.name = name;
DATA.mode = 'list';

j = 1; 
while j <= length(varargin)
    if strncmpi(varargin{j},'runallvpcs',7)
        DATA.mode = 'runallvpcs';
    end
    j = j+1;
end
%CheckDir([], DATA);
DATA.timerobj = timer('TimerFcn',{@CheckDir, DATA}, 'Period', 1.0,...
'ExecutionMode','FixedSpacing');
start(DATA.timerobj);

function CheckDir(a,b, DATA, varargin)

d = dir([DATA.name,'/*ClusterTimes.mat']);
for j = 1:length(d)
    name = [DATA.name '/' d(j).name];
    dname = strrep(name,'ClusterTimes','ClusterTimesDetails');
    dd = dir(dname);
%if ClusterTimes is newer that Details, needs recutting. But be sure that 
%Clustetimes is at least a minute old to leave time for peopleto make 2
%cuts
    if length(dd) ==1 && dd.datenum < d(j).datenum-0.001  && dd.datenum < now+0.001;
        if strcmp(DATA.mode,'list');
            fprintf('Need Cut for %s\n',name);
        elseif strcmp(DATA.mode,'runallvpcs')
            fprintf('Need Cut for %s\n',name);
            vname = strrep(name,'ClusterTimes','FullV');
            AllVPcs(vname,'tchan','all','savespikes','quantifyall','noninteractive');
        end
    end
end