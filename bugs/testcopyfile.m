function testcopyfile(laps)

sourcefilename = '/Volumes/bgc5/bgctestsourcefile.mat';

testdir = {'/Volumes/bgc5/smr/bgctargettestfile.mat' '/Volumes/bgc8/smr/bgctargettestfile.mat'  '/Volumes/bgc9/smr/bgctargettestfile.mat'};
if nargin == 0
    laps = 10;
end
if ischar(laps)
    testdir{1} = laps;
end


for j = 1:length(testdir)
    sometargetname = testdir{j};
    fprintf('Testing %s\n',testdir{j});
    delete(sometargetname); 
    for i = [1:laps]; 
        disp(i); 
        %kashk.r = randi(10, [10000000,1]); 
        copyfile(sourcefilename, sometargetname); 
        load(sometargetname); 
    end
end
