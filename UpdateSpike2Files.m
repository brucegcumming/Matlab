function UpdateSpike2Files(varargin)
%UpdateSpike2Files Copys s2s scripts from network to local drives
%run on Spike2 Experimental PCs

copyscript('/bgc/spike2/scripts/bgcs devel.s2s','C:/Spike2/bgcs.2s');
copyscript('/bgc/spike2/scripts/slave.s2s','C:/Spike2/slave.s2s');
copyscript('/bgc/spike2/scripts/SetConfig.s2s','C:/Spike2/SetConfig.2s');
copyscript('/bgc/spike2/scripts/MakeSpikeWaves.s2s','C:/Spike2/MakeSpikeWaves.2s');
copyscript('/bgc/spike2/scripts/MakeMaton.s2s','C:/Spike2/MakeMaton.2s');


function copyscript(src,tgt)

a = dir(src);
b = dir(tgt);
go = 1;
if isempty(a)
    cprintf('Cant find %s\n',src);
    return;
end
if ~isempty(a) && ~isempty(b) && b.datenum > a.datenum %target newer
    go = confirm(sprintf('%s (%s) is newer than %s (%s). Sure you want to copy?',tgt,b.date,src,a.date));
end
if ~isempty(a) && ~isempty(b) && b.datenum == a.datenum %target newer
    cprintf('blue','%s is already up to date (%s)\n',tgt,b.date);
    go = 0;
end
if go
    fprintf('Copying %s to %s\n',src,tgt);
    copyfile(src,tgt);
end
