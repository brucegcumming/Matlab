function [V, DATA] = ReadSpikeFiles(DATA, name)%[DATA, V] = ReadSpikeFiles(DATA, name) Read In Spike Files%if DATA is an integer, use that as probe number    V.name = name;    if isnumeric(DATA)        p = DATA;        clear DATA;        DATA.probe = p;        DATA.name = name;        DATA.clustersubdir = [];        DATA.exptno = GetExptNumber(DATA.name);    else        id = regexp(V.name,'Expt[0-9]*');        DATA.exptno = sscanf(V.name(id+4:end),'%d');    end    DATA.Expt.exptno = DATA.exptno;    s = AllV.SpkFileName(DATA);    V.Spikes = ReadSpikeFile(s,'allprobes');    V.exptno = DATA.exptno;    [monk, monkey, mdir] = GetMonkeyName(DATA.name);    if isempty(V.Spikes)        return;    end    if V.Spikes.Header.bysuffix == 0    if isfield(V.Spikes,'matfile')        V.matfile = V.Spikes.matfile;    else        [a,b] = fileparts(DATA.name);        V.matfile = [a '/' monk mdir '.mat'];    end    end