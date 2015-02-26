function [DATA, ispk] = APlotTrialSpikes(DATA, times, colors, nclusters, classify, varargin)

j = 1;
Trial = [];
ispk = [];
probe = DATA.probe;
lineoff = 0;
step = DATA.plot.SpikeMaxV;
step = DATA.vstep;
timemode = 0;
nprobes = 1;
plotxy = 1;
if DATA.plot.showartifacts
    maxcl = 9;
else
    maxcl = 8;
end
voff = NaN;
vscale = 1;
ip = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'lineoff',6)
        j = j+1;
        lineoff = varargin{j};
    elseif strncmpi(varargin{j},'Probes',6)
        j = j+1;
        probe = varargin{j};
    elseif strncmpi(varargin{j},'Probe',4)
        j = j+1;
        probe = varargin{j}(1);
        if length(varargin{j}) > 1
            ip = varargin{j}(2); % # in list
        end
        if length(varargin{j}) > 2
            nprobes = varargin{j}(3);
        end
    elseif strncmpi(varargin{j},'timemode',4)
        timemode = 1;
    elseif strncmpi(varargin{j},'noxy',4)
        plotxy = 0;
    elseif strncmpi(varargin{j},'spkid',5)
        j = j+1;
        ispk = varargin{j};
    elseif strncmpi(varargin{j},'Trial',4)
        j = j+1;
        Trial = varargin{j};
    elseif strncmpi(varargin{j},'vscale',4)
        j = j+1;
        vscale = varargin{j};
    elseif strncmpi(varargin{j},'voff',4)
        j = j+1;
        voff = varargin{j};
    end
    j = j+1;
end
if isfield(DATA,'AllSpikes')
    Spks = DATA.AllSpikes{probe(1)};
    if isfield(Spks,'pcs')
        PCs = Spks.pcs;
    else
        PCs = [];
    end
elseif DATA.probe == 100
    Spks = DATA.AllData.UstimV;
    PCs = Spks.codes;
else
    Spks = DATA.AllData.Spikes;
    if isempty(DATA.AllData.pcs)
        PCs = DATA.AllData.Spikes.codes;
    else
        PCs = DATA.AllData.pcs;
    end
end

if DATA.state.recut && size(Spks.codes,2) > 1
    ctype = 2;
else
    ctype = 1;
end

if ~isfield(Spks,'values')
    ispk = [];
    return;
end
splen = size(Spks.values,2).*size(Spks.values,3);
% try calculating energy for all spikes in Expt in one step.
if isempty(ispk)
    ispk = find(Spks.times > times(1) &...
        Spks.times < times(2) & Spks.codes(:,ctype) < maxcl);
end
if DATA.TriggerSign == -1
    id = find(Spks.values(ispk,11) < 0);
    ispk = ispk(id);
end
if isfield(DATA,'AllClusters') && DATA.bysuffix == 0 && DATA.state.somespikes == 0
    return;
end
if ispk
    if length(PCs) < max(ispk)
        PCs = Spks.codes;
    end
    if isfield(DATA,'sids')
        syncspk = find(ismember(ispk,DATA.sids{ip}));
        syncspklst = ismember(ispk,DATA.sids{ip});
        if length(probe) == 2
            Spks.values = (DATA.AllSpikes{probe(1)}.values(DATA.sids{1},:)+DATA.AllSpikes{probe(2)}.values(DATA.sids{2},:))/2;
            Spks.times =  DATA.AllSpikes{probe(1)}.times(DATA.sids{1});
            ispk = find(Spks.times > times(1) &...
                Spks.times < times(2));
        end
        if ismember(DATA.syncsign,[3 4 5])
            Spks.codes(ispk,3) = 0;
            Spks.codes(ispk(syncspk),3) = 1;
            ctype = 3;
        end
    end
    
    if DATA.probe < 100
        if isfield(DATA,'AllClusters')
            pid = GetProbe(DATA, DATA.currentexpt(1), DATA.probe);
            cx = DATA.AllClusters{DATA.currentexpt(1)}(pid).cx(ispk);
            cy = DATA.AllClusters{DATA.currentexpt(1)}(pid).cy(ispk);
        else
            [cx, DATA] = GetSpikeVals(DATA,ispk, Spks.values(ispk,:), Spks.dVdt(ispk,:),DATA.plot.clusterX, classify,PCs(ispk,:));
            DATA.Spikes.cx(ispk) = cx;
            %      [cy, DATA] = GetSpikeVals(DATA,ispk, SPKVARE, classify);
            [cy, DATA] = GetSpikeVals(DATA,ispk, Spks.values(ispk,:), Spks.dVdt(ispk,:),DATA.plot.clusterY, classify,PCs(ispk,:));
            DATA.Spikes.cy(ispk) = cy;
            dvdt = Spks.dVdt(ispk,:);
        end
        if length(ispk)
            DATA.currentspike = ispk(1);
        end
        adc = Spks.values(ispk,:);
        % recut == 2 means that clusters are not set here, but clusters have
        % been defined (previous list), so use those properites.
        if DATA.state.recut == 2 && DATA.probe == probe
            DATA = SetSpkCodes(DATA,ispk, DATA.probe,0);
            if isfield(DATA,'AllSpikes') %copy any changes into spks
                Spks.codes = DATA.AllSpikes{DATA.probe}.codes;
            else
                Spks.codes = DATA.AllData.Spikes.codes;
            end
        end
    else
        Spks.codes = DATA.AllData.UstimV.codes;
        adc = Spks.values(ispk,:);
        cx = zeros(size(ispk));
        cy = zeros(size(ispk));
    end
 adc(:,end+1) = NaN;
 xpts = 1:(splen+1);
 %method 2 MUCH faster in 2014b
method = 2;    
    for j = [1:nclusters+1]
        if method ==2
            cid = find(Spks.codes(ispk,ctype) == j-1);
            vs{j} = reshape(adc(cid,:)',[],1);
            xs{j} = reshape(repmat(xpts,1,length(cid)),1,[]);
        else
            vs{j} = [];
            xs{j} = [];
        end
    end
if method == 1    
    for spk = 1:length(ispk);
        if DATA.syncsign < 2 & isfield(DATA,'sids')
            j = syncspklst(spk)+1;
        else
            j = Spks.codes(ispk(spk), ctype)+1;
        end
        if DATA.plot.dvdt ==2
            vs{j} = [vs{j} dvdt(spk,:) NaN];
            xs{j} = [xs{j} adc(spk,1:splen-1) NaN];
        elseif DATA.plot.dvdt
            vs{j} = [vs{j} dvdt(spk,:) NaN];
            xs{j} = [xs{j} [1:splen-1] NaN];
        elseif timemode
            vs{j} = [vs{j} adc(spk,:) NaN];
            xs{j} = [xs{j} [1:splen]+(Spks.times(ispk(spk))-times(1)).*timemode NaN];
        elseif DATA.plot.nodc
            vs{j} = [vs{j} adc(spk,:)-mean(adc(spk,:)) NaN];
            xs{j} = [xs{j} [1:splen] NaN];
        elseif DATA.plot.voltxy == 3
            vs{j} = [vs{j} adc(spk,65:96) NaN];
            xs{j} = [xs{j} adc(spk,97:128) NaN];
        elseif DATA.plot.voltxy == 2
            vs{j} = [vs{j} adc(spk,33:64) NaN];
            xs{j} = [xs{j} adc(spk,65:96) NaN];
        elseif DATA.plot.voltxy == 1
            vs{j} = [vs{j} adc(spk,1:32) NaN];
            xs{j} = [xs{j} adc(spk,33:64) NaN];
        elseif DATA.plot.voltxy == 4
            vs{j} = [vs{j} adc(spk,65:96)-adc(spk,33:64) NaN];
            xs{j} = [xs{j} [1:splen/4] NaN];
        else
            vs{j} = [vs{j} adc(spk,:) NaN];
            xs{j} = [xs{j} xpts NaN];
        end
    end
end
    if length(DATA.plot.voffsets) >= ip
        for j = [1:nclusters+1]
            vs{j} = vs{j}+DATA.plot.voffsets(ip);
        end
        DATA.plot.currentoffset(ip) = DATA.plot.voffsets(ip);
    else
        DATA.plot.currentoffset(ip) = (ip-(nprobes+1)/2).*step;
        for j = [1:nclusters+1]
            vs{j} = vs{j}+(ip-(nprobes+1)/2).*step;
        end
    end
    
    if timemode && isfield(DATA,'timefig')
        set(0,'CurrentFigure',DATA.timefig);
        text(0,(ip-(nprobes+1)/2).*step,sprintf('%d',probe));
        vh = DATA.tvh;
    else
        set(0,'CurrentFigure',DATA.svfig);
        vh = DATA.svh;
    end
    nc = min([nclusters+1 length(DATA.svh)]);
    
    
    for j = 1:nc
        k = j+lineoff;
        if k > length(vh)
            vh(k) = line('Xdata' , xs{j}, 'Ydata', vs{j});
        elseif ~isempty(xs{j})
            if ishandle(vh(k))
                set(vh(k),'Xdata' , xs{j}, 'Ydata', vs{j});
            else
                vh(k) = line('Xdata' , xs{j}, 'Ydata', vs{j});
            end
        elseif double(vh(k)) > 0 && ishandle(vh(k)) %no spikes, wipe clean
            set(vh(k),'Xdata' , 0, 'Ydata', 0);
        end
    end
    if ~isempty(Trial)
        if length(DATA.probelist) > 2
            xc = sprintf('P%d',DATA.probe);
        else
            xc = [];
        end
        xc = [xc cmb.ExtraLabels(Trial)];
        nspk = sum(Spks.codes(ispk,ctype) == DATA.currentcluster);
        title(sprintf('Trial %d (id%d %.2f - %.2f) ed%.3f%s %d/%d spks',abs(Trial.Trial),...
            Trial.id,Trial.Start(1)./10000,Trial.End(end)./10000,Trial.ed,xc,nspk,length(ispk)));
    end
    
    if probe == DATA.probe && plotxy
        drawnow;
        set(0,'CurrentFigure',DATA.xyfig);
        for j = 0:nclusters
            sp = find(Spks.codes(ispk, ctype) == j);
            plot(cx(sp),cy(sp),...
                '.','color',DATA.spkcolor{j+1},'markersize',DATA.ptsize);
            hold on; %% need this when called from PlotOneTrial
        end
    end
else
    title(sprintf('Trial %d (id%d %.2f - %.2f) ed%.3f',abs(Trial.Trial),...
        Trial.id,Trial.Start(1)./10000,Trial.End(end)./10000,Trial.ed));
    
end
DATA.minplottime = 0.00;
if DATA.minplottime > 0.001
    while toc < DATA.minplottime
    end
end


