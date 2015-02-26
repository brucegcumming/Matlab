function FixSpkBlk(path, timefix)
% FixSpkBlk(namestr, fixtime)
% fixes time mismatch for spkblk files caused by  trigger errors 
% at runtime in SPike2.
% namestr is passed to dir to get a list of names
% 
% make the search as specific as possible
% FixSpkBlk('/Volumes/bgc9/smr/lem/M270/lemM270.1.spkblk*',fixtime)
%
%if fixtime is 0, re-writes all the files with no fixtime

d = mydir(path);
    load(d(j).name);
    if timefix ~= 0
    save(d(j).name,'Ch*','timefix');
    else
    save(d(j).name,'Ch*'); %save without fixtime to undo.
    end
        
end