function testsave(laps)

if nargin == 0
    laps = 10;
end

somefilename = '/Volumes/bgc5/smr/bgctestfile.mat'; %mat4
delete(somefilename); for i = [1:laps]; disp(i); kashk.r = randi(10, [10000000,1]); save(somefilename, 'kashk'); load(somefilename); end


%somefilename = '/Volumes/bgc9/smr/bgctestfile.mat'; %ds1
somefilename = '/Volumes/bgc9/testnewdir2/bgctestfile.mat'; %ds1
delete(somefilename); for i = [1:laps]; disp(i); kashk.r = randi(10, [10000000,1]); save(somefilename, 'kashk'); load(somefilename); end

