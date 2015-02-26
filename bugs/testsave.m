function testsave(laps, varargin)

delay = 0;
deleteprev = 0;
host = [];
testdir = {'/Volumes/bgcdiv/data/test/bgctestfile.mat' '/Volumes/bgcpriv/tmp/bgctestfile.mat'};

delays(1) = 0.0;
j = 1;
while j <= length(varargin)
    if sum(strcmp(varargin{j},{'mat4','priv','tmp'}))
        host = varargin{j};
    elseif strcmp(varargin{j},'delete')
        deleteprev = 1;
    elseif strcmp(varargin{j},'move')
        deleteprev = 1;
    elseif strcmp(varargin{j},'delay')
        j = j+1;
        delays = varargin{j};
    end
    j = j+1;
end


if ~isempty(host)
    if strcmp(host,'mat4')
        testdir = {'/Volumes/bgc7/bgc8/smr/testdir/bgctestfile.mat' '/Volumes/bgc7/smr/bgctestfile.mat' '/Volumes/bgc5/smr/bgctestfile.mat'};
    elseif strcmp(host,'priv')
        testdir = {'/Volumes/bgcpriv/tmp/bgctestfile.mat'};
    elseif strcmp(host,'tmp')
        testdir = {'/Volumes/tmp/bgctestfile.mat'};
    end
end

if nargin == 0
    laps = 10;
end
if ischar(laps)
    testdir{1} = laps;
    laps = 10;
elseif iscellstr(laps)
    testdir = laps;
    laps = 10;
end

for j = 1:length(testdir)
somefilename = testdir{j};
fprintf('Testing %s\n',testdir{j});
backupname = strrep(somefilename,'.mat','bak.mat');
for i = [1:laps]; 
%either method of deleting the file first seems to prevent the problems
if exist(somefilename) && deleteprev ==1
        delete(somefilename);
    elseif exist(somefilename) && deleteprev ==2 %use movefile
        movefile(somefilename,backupname);
    end
    kashk.r = randi(10, [10000000,1]); 
    fprintf('Saving %d: %s at %s\n',i,somefilename,datestr(now));
    %-v7.3  produces error at write time.
%-v6  gets rid of error
    save(somefilename, 'kashk');
    pause(delays(1));
    clear('kashk');
    fprintf('Loading %s at %s\n',somefilename,datestr(now));
    load(somefilename); 
     pause(delays(1));
end
delete(somefilename); 
end

