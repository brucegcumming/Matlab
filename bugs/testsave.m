function testsave(laps)

testdir = {'/Volumes/bgc5/smr/bgctestfile.mat' '/Volumes/bgc9/smr/bgctestfile.mat'};
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
delete(somefilename); for i = [1:laps]; kashk.r = randi(10, [10000000,1]); 
    fprintf('Saving %d: %s at %s\n',i,somefilename,datestr(now));
    save(somefilename, 'kashk'); 
    fprintf('Loading %s at %s\n',somefilename,datestr(now));
    load(somefilename); end
end

