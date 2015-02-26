function RestoreFile(name, varargin)
%RestoreFile(name, ...    Restore Backups made with BackupFile
%Finds most recent file in backup directory matching name and copies
%back to the original target.

dpath = fileparts(name);
bpath = strrep(name,dpath,[dpath '/backup/']);
bpath = strrep(bpath,'.mat','*');
d = dir(bpath);
[a,b] = sort([d.datenum],'descend');

broot = fileparts(bpath);
for j = 1:length(b)
    yn = questdlg(sprintf('Copy %s to %s', d(b(j)).name,name), 'File Restore', 'Yes', 'Next', 'Cancel','Yes');
    if strcmpi(yn,'yes')
        copyfile([broot '/' d(b(j)).name],name);
        return;
    elseif strcmpi(yn,'Cancel')
        return;
    end
end