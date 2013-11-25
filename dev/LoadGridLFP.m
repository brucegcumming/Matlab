function Expt = LoadGridLFP(Expt, varargin)
preperiod = 1000;
postperiod = 1000;
verbose = 0;
if isfield(Expt.Stimvals,'po')
    postperiod = Expt.Stimvals.po .* 10000;
end
if isfield(Expt.Stimvals,'pr')
    preperiod = Expt.Stimvals.pr .* 10000;
end

zfill = 0;
needft = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'postperiod',5)
        j = j+1;
        postperiod = varargin{j} .* 10000;
    elseif strncmpi(varargin{j},'ft',2)
        needft = 1;
    elseif strncmpi(varargin{j},'fillzero',6)
        zfill = 1;
    elseif strncmpi(varargin{j},'preperiod',6)
        j = j+1;
        postperiod = varargin{j} .* 10000;
    elseif strncmpi(varargin{j},'verbose',6)
        verbose = 1;
    end
    j = j+1;
end

Expt.Header.preperiod = preperiod;
Expt.Header.postperiod = postperiod;

loaddir = fileparts(Expt.Header.loadname);
eid = [];
if isfield(Expt.Header,'Combineids')
    eid = Expt.Header.Combineids;
elseif isfield(Expt.Header,'idrange')
    idx = BuildGridIndex(fileparts(Expt.Header.loadname),[]);
    eid = find(idx.exptids == Expt.Header.idrange(end));
end
gap = [];
if ~isempty(eid)
    for j = 1:length(Expt.Trials)
        Expt.Trials(j).lfplen = 0;
    end
    for j = 1:length(eid)
        name = sprintf('%s/Expt%d.lfp.mat',loaddir,eid(j));
        if exist(name,'file')
        ts = now;
        if verbose
            fprintf('Loading %s....',name);
        end
        load(name);
        LFP.rawlfp = double(LFP.rawlfp) ./ LFP.intscale(1);
        LFP.t = LFP.t .* 10000;
        Expt.Header.LFPsamplerate = LFP.samper;
        Expt.Header.CRsamplerate = LFP.samper/10;
        Expt.Header.LFPsigma = LFP.sd;
        Expt.Header.LFPdecimate = LFP.decimate;
        for j = 1:length(Expt.Trials)
            id = find(LFP.t > Expt.Trials(j).Start(1)-preperiod & ...
                LFP.t < Expt.Trials(j).End(end)+postperiod);
            if ~isempty(id)
                offset = round((LFP.t(id(1))- (Expt.Trials(j).Start(1)-preperiod))./(LFP.samper .*10000));
                id = id-offset;
                gap(j) = LFP.t(id(1))-LFP.t(id(1)-1);
                Expt.Trials(j).lfpstart = LFP.t(id(1))-Expt.Trials(j).Start(1);
                Expt.Trials(j).lfptime = LFP.t(id(1));
                Expt.Trials(j).lfplen = length(id);
                Expt.Trials(j).LFP = LFP.rawlfp(id,:);
                if needft
                    Expt.Trials(j).FTlfp = abs(fft(Expt.Trials(j).LFP));
                end
            end
        end
        if verbose
            fprintf('%.1f\n',mytoc(ts));
        end
        else
            mycprintf('red','Cant Read %s\n',name);
        end
    end
end
if ~isempty(gap)
    x = std(gap);
end