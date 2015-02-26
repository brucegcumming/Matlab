function Spikes = GetProbeFiles(DATA, probe, subprobe, varargin)
%GetProbeFiles loads up spikes that are spli across files.
% given a time range (in sec, not tics) only loads files needed for that
% range
trange = [];
dvfile = [];
nodv = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j}, 'trange',3)
        j = j+1;
        trange = varargin{j};
    elseif strncmpi(varargin{j}, 'nodv',3)
        nodv = 1;
    end
    j = j+1;
end
id = find([DATA.probes.probe] == probe);
[dp,pref] = fileparts(DATA.datafilename);
Spk.times = [];
Spk.values = [];
Spk.codes = [];
[a,sid] = sort([DATA.probes(id).first]);
sid = id(sid);
if length(trange) > 1
    %fid is files that are past the point we need
    fid = find([DATA.probes(sid).first] > trange(2));
    if isempty(fid)
        fid = length(sid)+1;
    end
    lid = find([DATA.probes(sid).last] < trange(1));
    if isempty(lid) %need first file
        sid = sid(1:fid(1)-1);
    elseif fid(1) == lid(end)+1
        sid = sid(fid(1));
    else
        sid = sid(lid(end)+1:fid(1)-1);
    end
    fprintf('Spike file %d\n',sid);
end
if DATA.state.savedvdt
    dvfile = [dp '/Spikes/' strrep(DATA.probes(sid(1)).filename,'t0','dvdt')];
end
for j = 1:length(sid)
    filename = [dp '/Spikes/' DATA.probes(sid(j)).filename];
    if ~exist(filename,'file') & isfield(DATA.probes,'pathname');
        filename = [DATA.probes(sid(j)).pathname '/' DATA.probes(sid(j)).filename];
    end
    if exist(filename,'file')
        a = load(filename);
        chname = DATA.probes(sid(j)).var;
        if isempty(Spk.times)
            Spk = a.(chname);
            ns = length(Spk.times);
        else
            Spk.times = [Spk.times; a.(chname).times];
            Spk.codes = [Spk.codes; a.(chname).codes];
            Spk.values = [Spk.values; a.(chname).values];
            if size(a.(chname).times,1) > size(a.(chname).values,1)
                fprintf('Some Missing Spike values in %s\n',filename);
                if DATA.logfid
                    fprintf(DATA.logfid,'Some Missing Spike values in %s\n',filename);
                end
            end
        end
    else
        fprintf('No file %s\n',filename);
    end
end

if ~isempty(Spk.values)
    if size(Spk.values,1) < length(Spk.times)
        fprintf('Some Missing Spike values in %s\n',pref);
    end
    
    if size(Spk.values,3) > 1 & subprobe > 0
        Spk.values = Spk.values(:,:,subprobe);
    end
    if nodv == 0
        Spikes = CleanSpikes(Spk,'bufl',10000,'dvfile',dvfile);
    else
        Spikes = Spk;
    end
    Spikes.times = Spikes.times .* 10000;
    if isfield(DATA.probes, 'firsti')
        Spikes.firstspki = DATA.probes(sid(1)).firsti;
    end
    cfile = cmb.ClusterFile(DATA,'probe',probe);
    if exist(cfile,'file')
        load(cfile);
        lastspk = DATA.probes(sid(end)).firsti+DATA.probes(sid(end)).nspk-1;
        if length(clid) < lastspk
            lastspk = length(clid);
            nspk = 1+lastspk-Spikes.firstspki;
            Spikes.codes(1:nspk,2) = clid(Spikes.firstspki:lastspk);
        elseif size(Spikes.codes,1) >  lastspk
            Spikes.codes(Spikes.firstspki:lastspk,2) = clid(Spikes.firstspki:lastspk);
        else
            Spikes.codes(:,2) = clid(Spikes.firstspki:lastspk);
        end
        Spikes.codes = Spikes.codes(:,1:2);
    end
else
    fprintf('No Spike values %s\n',filename);
    Spikes = [];
end



