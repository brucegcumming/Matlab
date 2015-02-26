function NEV = ns2nev(name, varargin)

nsmp = 32;
savenev = 1; %save if new
pre = 8;
th = -100;
NEV = [];
rebuild = 0;

[spkdir, nsname] = fileparts(name);
spkdir = [spkdir '/GridSpikes'];

nevfile = strrep(name,'.ns5','.nev');
if exist(nevfile,'file'); %% get Events, spike trigger levels.
    NEV = openNEV('read',nevfile);
end

if isempty(NEV.Data.Spikes.Waveform) || rebuild
    NEV.Data.Spikes = [];
nsx = openNSx('read',name);
for e = 1:size(nsx.Data,1)
    v = nsx.Data(e,:);
    th = std(v) * -3;
    t = find(v(nsmp:end-nsmp) >= th & v(nsmp-1:end-nsmp-1) < th)' + nsmp-1;
    ts = repmat(t,1,nsmp) + repmat([1:nsmp]-pre,length(t),1);
    ts(find(ts <1)) = 1;
    if e ==1
        NEV.Data.Spikes.Waveform = v(ts);
        NEV.Data.Spikes.Electrode = ones(size(t));
        NEV.Data.Spikes.TimeStamp = t;
    else
        NEV.Data.Spikes.Waveform = cat(1,NEV.Data.Spikes.Waveform,v(ts));
        NEV.Data.Spikes.Electrode = cat(1, NEV.Data.Spikes.Electrode, ones(size(t)) .* e);
        NEV.Data.Spikes.TimeStamp = cat(1, NEV.Data.Spikes.TimeStamp, t);
    end
end
NEV.Data.Spikes.Unit = zeros(size(NEV.Data.Spikes.Electrode));

matfile = strrep(name,'.ns5','.mat'); %if we had to re-read the nsx file, we need to re-write the mat file
if savenev
    save(matfile,'NEV');
end

elseif exist(spkdir,'dir')
    d = dir(spkdir);
    np = 0;
    for j = 1:length(d)
        if strncmp(nsname, d(j).name, length(nsname))
            a = load([spkdir '/' d(j).name]);
            np = np+1;
            if isfield(a.nev,'probeid')
                p = a.nev.probeid;
            else
                p = sscanf(d(j).name(2+length(nsname):end),'%d');
            end
            t = a.nev.Data.Spikes.TimeStamp;
            if np == 1
                NEV.Data.Spikes.Waveform = a.nev.Data.Spikes.Waveform;
                NEV.Data.Spikes.Electrode = ones(size(t)).*p;
                NEV.Data.Spikes.TimeStamp = t;
            else
                NEV.Data.Spikes.Waveform = cat(1,NEV.Data.Spikes.Waveform,a.nev.Data.Spikes.Waveform);
                NEV.Data.Spikes.Electrode = cat(1, NEV.Data.Spikes.Electrode, ones(size(t)) .* p);
                NEV.Data.Spikes.TimeStamp = cat(1, NEV.Data.Spikes.TimeStamp, t);
            end
        end
    end
    NEV.Data.Spikes.Unit = ones(size(NEV.Data.Spikes.Electrode));
end