function png2pgm(dname, varargin)

overwrite = 0;
outdir = [];
saveim = 1;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'outdir',6)
        j = j+1;
        outdir = varargin{j};
    elseif strncmpi(varargin{j},'overwrite',6)
        overwrite = 1;
    end
    j = j+1;
end

d = dir([dname '/IM*.png']);
for j = 1:length(d)
    name = [dname '/' d(j).name];
    if isempty(outdir)
        outname = [dname '/' d(j).name];
    else
        outname = [outdir '/' d(j).name];
    end
        
    task = 'Writing';
    if regexp(d(j).name,'\.png$')
        outname = strrep(outname,'png','pgm');
        if exist(outname) && ~overwrite
            fprintf('%s already exists\n',outname);
        task = 'OverWriting';
        elseif saveim  
            fprintf('Reading %s\n',name);
            im = imread(name);
            imwrite(im,outname,'pgm');
        end
    end
end