function BackupFile(name, varargin)
% BackupFile(name)
% Takes a full  path name 'path/file', and moves the named file to 'path/backup/filemmddyy.HHMMSS' (date.time of backup  in name)
% BackupFile(name,'copy') copies the file rather than moving. Slower, but safer if you are not about to overwrite
% BackupFile(,...'showtime', uses mmddyy.HHMMSS, includeing current time and date
% in the backupname. This is the default.
% BackupFile(,...'notime', uses mmddyy, includeing time and date in the backupname
% BackupFile(,...'filetime', uses the files time/date in the name
% creates backup directory if absent
% if name refers to a nonexistent file, no action is taken
suffix = [];
docopy = 0;
verbose = 0;
usefiletime = 0;
datefmt = 'mmddyy.HHMMSS';
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'suffix',4)
        j = j+1;
        suffix = varargin{j};
    elseif strncmpi(varargin{j},'notime',6)
        datefmt = 'mmddyy';
    elseif strncmpi(varargin{j},'filetime',6)
        usefiletime = 1;
    elseif strncmpi(varargin{j},'showtime',6)
        datefmt = 'mmddyy.HHMMSS';
    elseif strncmpi(varargin{j},'copy',4)
        docopy = 1;
    elseif strncmpi(varargin{j},'print',4)
        verbose = 1;;
    end
    j = j+1;
end
[a,b,c] = fileparts(name);
if isempty(a)
    backdir = './backup';
else
    backdir = [a '/backup'];
end
if ~exist(name,'file')
    return;
end
if ~exist(backdir,'dir')
    mkdir(backdir);
end
if exist(backdir,'dir')
    if usefiletime
        d = dir(name);
        dstr = datestr(d.datenum,datefmt);
    else
        dstr = datestr(now,datefmt);
    end
    backfile = [backdir '/' b dstr suffix c];
    ts = now;
    if docopy
        str = 'Copying';
        copyfile(name, backfile);
    else
        try
        movefile(name, backfile);
        catch ME
            CheckExceptions(ME);
            cprintf('red','Error Trying to Move %s to %s\n',name,backfile);
        end
        str = 'Moving';
    end
    if verbose
        fprintf('%s %s to %s (%.2f sec)\n',str,name,backfile,mytoc(ts));
    end
else
    cprintf('errors','Can''t Backup to %s - Folder nonexistent\n',backdir);
end
