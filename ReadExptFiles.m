function [AllExpts, xcorrs] = ReadExptFiles(path, prefix)

AllExpts = {};
xcorrs = [];

d = dir(path);
nx = 0;
for j = 1:length(d)
    if regexp(d(j).name,[prefix '.*cell[0-9]*'])
        name = [path '/' d(j).name];
        load(name);
        nx = nx+1;
        AllExpts{nx} = Expt;
    end
    if regexp(d(j).name,[prefix '.*xcorrs'])
        name = [path '/' d(j).name];
        load(name);
    end
end