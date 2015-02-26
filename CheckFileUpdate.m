function CheckFileUpdate(src, tgt, chkmode)
%CheckFileUpdate(src, tgt, chkmode) Copy new files
%checks if tgt is older than src, and pops up a window
%to allow copy
     if nargin < 3
        chkmode = 'change';
     end
    a = dir(src);
    b = dir(tgt);
    if isempty(a)
        fprintf('Missing %s - cant check update\n',src);
        return;
    end
    if isempty(b) %target does not exist - just copy
            try  %This often produces permission error
                copyfile(src,tgt);
            catch ME
                cprintf('errors',ME.message);
                fprintf('Error copying %s\n',tgt);
            end
    end
    if strncmp(chkmode,'new',3)         
        return;
    end
    if ~isempty(a) && ~isempty(b) && a.datenum > b.datenum
        yn = questdlg(sprintf('%s is newer. Copy to %s?',src,tgt),'Update Check','Yes','No','Yes');
        if strcmp(yn,'Yes')
            try  %This will produce and error becuase verg.m is in use. But the copy succeeds
                [a,b,c] = copyfile(src,tgt);
            catch ME
                cprintf('errors',ME.message);
                fprintf('possible error copying %s\n',tgt);
            end
        end
    end
 