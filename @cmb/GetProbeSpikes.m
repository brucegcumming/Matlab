function Spikes = GetProbeSpikes(All, filename, varname, probe)
if length(probe) > 1
    subprobe = probe(2);
else
    subprobe = 1;
end
gotspk = 0;
ts = now;
if isfield(All,'toplevel')
    Spk = getappdata(All.toplevel,'Spike2File');
    if isfield(Spk,'filename') && strcmp(Spk.filename,filename) 
        gotspk = 1;
    end
end
if gotspk
    
elseif exist(filename,'file')
    Spk = load(filename);
    Spk.filename = filename;
    if isfield(All,'toplevel')
        setappdata(All.toplevel,'Spike2File',Spk);
    end
    if isfield(Spk,'Spikes')
        Spikes = Spk.Spikes;
    end
else
    fprintf('No file %s\n',filename);
    Spikes = [];
    return;
end
loaddur = mytoc(ts);
if isempty(varname)
    Spikes.times = Spikes.times .* 10000;
    if size(Spikes.codes,2) == 1
        Spikes.codes(:,2) = Spikes.codes(:,1);
    end
    Spikes.loaddur = loaddur;
elseif isfield(Spk,varname)
    %        Spikes = eval(varname);
    Spikes = CleanSpikes(Spk.(varname), 'bufl',10000);
    Spikes.times = Spikes.times .* 10000;
    Spikes.loaddur = loaddur;
else
    fprintf('No data for %s in %s\n',varname,filename);
    Spikes = [];
    return;
end
if isinteger(Spikes.values) && isfield(Spikes,'maxv')
    if isfield(Spikes,'maxint')
        Spikes.values = double(Spikes.values) .* Spikes.maxv/Spikes.maxint;
    else
        Spikes.values = double(Spikes.values) .* Spikes.maxv/32000;
        Spikes.maxint = 32000;
    end
end
ispk = find(Spikes.codes(:,1) > 0); %exclude codes == 0 for some cases?

if size(Spikes.values,3) > 1 & subprobe > 0
    Spikes.values = Spikes.values(:,:,subprobe);
end
if isempty(probe)
    Spikes.probe = 0;
else
    Spikes.probe = probe(1);
end

