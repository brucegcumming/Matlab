function [im, details] = ReadPGM(name, varargin)
%[im, details] = ReadPGM(name... read pgm image and comments

details = [];
im = imread(name);
fid = fopen(name,'r');
if fid > 0
    line = fgets(fid);
    ok = 1;
    while ok
        line = fgets(fid);
        if line(1) == '#'
            details.comments{ok} = deblank(line(2:end));
            ok = ok+1;
        else
            ok = 0;
        end
    end
    fclose(fid);
end
