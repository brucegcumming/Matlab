function SpkOnline(path, varargin)

d = dir(path)
for j = 1:length(d)
    if strncmp(d(j).name,'Expt',4)
        fstr