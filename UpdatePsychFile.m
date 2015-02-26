function UpdatePsychFile(name, varargin)
%UpdatePsychFile(name, copies local pysch file to network

ld = dir(name);
if isempty(ld)
    fprintf('No file %s to copy Psych from\n',name);
    return;
end
netname = strrep(name, '/local','/b')
d = dir(netname);
if isempty(d) || d.datenum < ld.datenum
    BackupFile(netname,'copy');
    success = copyfile(name, netname);
    if success
        fprintf('Copied %s to %s\n',name,netname);
    else
        cprintf('red','Error Copying %s to %s\n',name,netname);
    end
else        
    fprintf('%s is up to date\n',netname);
end
