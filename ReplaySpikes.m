function ReplaySpikes(lstname)

ix = regexp(lstname,'[0-9]\.[0-9]\.');
if ~isempty(ix)
    uffname = lstname(1:ix(1)+2);
end

type = computer;
if strmatch(computer,'PCWIN')
    olddir = pwd;
    cd('C:');
    cmd = sprintf('C:/Cygwin/bin/bash -login -c "/home/bgc/bin/Cygwin/replay -ms %s -f %s"',uffname,lstname);
    fprintf('Running %s\n',cmd);
    system(cmd);
    cd(olddir)
else
    cmd = sprintf('replay -ms %s -f %s',uffname,lstname);
    fprintf('Running %s\n',cmd);
    system(cmd);
end



%
%Under Windows need:
%
%bash -login -c "/home/bgc/bin/Cygwin/replay -ms lstname"
%(and need exceed running...)
