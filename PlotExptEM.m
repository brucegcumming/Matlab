function [mev, details] = PlotExptEM(Expt, varargin)
%
% PlotExptEM plots eye movement data (and loads the data into Expt if
% needed)
%
% PlotExptEM(Expt, 'choice') returns mean eye positions by choice
% PlotExptEM(Expt, 'choicediff') returns difference in mean conjugate eye
% positions by choice (DOWN/RIGHT - UP/LEFT)
% PlotExptEM(Expt, 'choicediff',pidx,nidx) gives ids for trials to calculate difference for
% (difference is nidx - pidx)
% In both cases, there is no adjustment for the neurons prefereed.   
%
% PlotExptEM(Expt, 'choicediff', 'left') just uses left eye data (later rufus files)
% PlotExptEM(Expt, 'choicediff', 'sumsac') looks at the vector sum of saccades in each trial, and
%            finds a subset of trials for which thie mean vector is the same for each choice
% PlotExptEM(Expt, 'choicediff', 'sumsac','stattime',[a b c d]) defines a time window for calculations
%             for saccade summing, all saccades after a and before d are included
% PlotExptEM(Expt, 'choicediff', 'sumsac','cptime',[a b]) only use spikes between a and b to calculate CP

%load('C:/bgc/data/psych/EmExpt')
%PlotExptEM(Expt,'choicediff','left','cmppos','stattime',[00 1000 19000 20000]); set(gca,'ylim',[-0.5 0.5])
%is use just last 100ms and make equal, get pretty flat lines, with a bit
%of null shift at start. But if select for small shifts whithin trial (compare start and end), get
%get large intial and final shifts, with zero in the middle!!
%PlotExptEM(Expt,'choicediff','left','cmpsac','stattime',[00 1000 19000 20000]); set(gca,'ylim',[-0.5 0.5])

trid = 1:length(Expt.Trials);
showplot = 1;
startzero = 0;
details.n = 0;
splitchoice = 0;
plotdiff = 0;
duration = 20000;
emskip = 0;
usechan = [1:4];
varargon = {};
em.maxdrift = 0;
em.stattime = [1500 6000 15500 20000];
em.hchan = [1 2];
em.vchan = [3 4];
em.cptime = [];
labelb = 'Saccades';
labela = 'MeanEyePos';

j = 1;
while j <= nargin-1
    str = varargin{j};
    if strncmpi(str,'choice',3)
        splitchoice = 1;
        if length(varargin) > j+1 & isnumeric(varargin{j+1})
            j = j+1;
            pid = varargin{j};
            j = j+1;
            nid = varargin{j};
        else
            pid = find([Expt.Trials.RespDir] < 0);  %UP/LEFT
            nid = find([Expt.Trials.RespDir] > 0); %DOWN/RIGHT
        end
        if strfind(str,'diff')
            plotdiff = 1;
        end
    elseif strncmpi(str,'allverg',6)
        varargon = {varargon{:} varargin{j}};
    elseif strncmpi(str,'cptime',6) %define time intervals for start/end periods
        j = j+1;
        em.cptime = varargin{j};
    elseif strncmpi(str,'duration',3)
        j = j+1;
        duration = varargin{j};
    elseif strncmpi(str,'emskip',3)
        j = j+1;
        emskip = varargin{j};
    elseif strncmpi(str,'left',4)
        usechan = [1 3];
        varargon = {varargon{:} varargin{j}};
        em.hchan = 1;
        em.vchan = 3;
    elseif strncmpi(str,'mean',4)
        varargon = {varargon{:} varargin{j}};
    elseif strncmpi(str,'sumsac',6)
        varargon = {varargon{:} varargin{j}};
    elseif sum(strncmpi(str,{'cmpsac' 'cmppos' 'startpos'},6))
        varargon = {varargon{:} varargin{j}};
    elseif strncmpi(str,'stattime',6) %define time intervals for start/end periods
        j = j+1;
        em.stattime = varargin{j};
    elseif strncmpi(str,'signal',3)
        j = j+1;
        signals = varargin{j};
        allid = [];
        if ~isfield(Expt.Trials,'signal')
            Expt = FillTrials(Expt,'signal');
        end
        for k = 1:length(signals)
            id = find([Expt.Trials.signal] == signals(k));
            allid = [allid id];
        end
        if ~isempty(allid)
        Expt.Trials = Expt.Trials(allid);
        end
        trid = 1:length(Expt.Trials);
    elseif strncmpi(str,'trials',3)
        j = j+1;
        trid = varargin{j};
    elseif strncmpi(str,'verg',4)
        varargon = {varargon{:} varargin{j}};
    elseif strncmpi(str,'zero',3)
        startzero = 1;
    end
    j = j+1;
end

if ~isfield(Expt.Trials,'EyeData')
    Expt = LoadEmData(Expt);
end
sample_rate = 1/(10 * Expt.Header.CRsamplerate);
eyevals = cat(3,Expt.Trials(trid).EyeData);
mev = mean(eyevals,3);
sds = std(mev);
if sds(2) < 1e-10 & sds(4) < 1e-10
    details.eyes = [1 0 1 0];
else
    details.eyes = [1 1 1 1];
end

times = ([1:size(mev,1)] .* 0.1./Expt.Header.CRsamplerate); %in ms
if isfield(Expt.Header,'emtimes')
times = times + Expt.Header.emtimes(1)./10;
times = Expt.Header.emtimes./10;
end
details.times = times;

if splitchoice & isfield(Expt.Trials,'Saccades')
 pm = PlotExptEM(Expt,'Trials',pid,varargon{:});
 psacs = [Expt.Trials(pid).Saccades];
 nsacs = [Expt.Trials(nid).Saccades];
 if isempty(psacs)
     details.prefsac = 0;
 else
 id = find([psacs.start] < duration & [psacs.start] > emskip);
 [x,y] = pol2cart([psacs(id).dir], [psacs(id).size]);
 details.prefsac = mean(x) + i * mean(y);
 details.prefsacs = x + i * y;
 details.prefsact = [psacs(id).start];
 end
 if isempty(nsacs)
     details.nullsac = 0;
 else
 id = find([nsacs.start] < duration & [nsacs.start] > emskip);
 [x,y] = pol2cart([nsacs(id).dir], [nsacs(id).size]);
 details.nullsac = mean(x) + i * mean(y);
 details.nullsacs = x + i * y;
 details.nullsact = [nsacs(id).start];
 end
 details.ntrials = [length(pid) length(nid)];
 nm = PlotExptEM(Expt,'Trials',nid,varargon{:});
 %remove from each choice the overall starting mean ACROSS all choices.
   sid = find(times > 100 & times < 300);
    for j = 1:size(mev,2)
        pm(:,j) = pm(:,j) - mean(mev(sid,j));
        nm(:,j) = nm(:,j) - mean(mev(sid,j));
    end
    mev = cat(3,pm, nm);
    if plotdiff
        mev = diff(mev,1,3);
        if details.eyes(2) == 0 %r eye bad
            emv(:,1) = mev(:,1);
            emv(:,2) = mev(:,3);
        else
            emv(:,1) = mean(mev(:,[1 2]),2);
            emv(:,2) = mean(mev(:,[3 4]),2);
        end
        mev = emv;
    end
else 
    pm = [];
end


if startzero
    for j = 1:size(mev,2)
        mev(:,j) = mev(:,j) - mean(mev(1:10,j));
    end
end

if showplot
    GetFigure(labela);
    hold off;
    if splitchoice & length(pm)
        if plotdiff
            plot(times,pm(:,usechan)-nm(:,usechan));
        else
            plot(times,pm(:,usechan));
                hold on;
            plot(times,nm(:,usechan),':');
        end
    else
        plot(times(1:size(mev,1)),mev(:,usechan));
        labels = {'LH','RH','LV','RV'};
        legend(labels{usechan});
    end
    if isfield(details,'prefsacs')
        GetFigure(labelb);
        sacres = PlotSaccByChoice(Expt.Trials, Expt.Header, pid,nid, em, varargon{:});
        details.sacres = sacres;
        if isfield(sacres,'adjpid')
        GetFigure(labela);
        mev = mean(eyevals(:,:,[sacres.adjpid sacres.adjnid]),3);
        pm = mean(eyevals(:,:,sacres.adjpid),3);
        nm = mean(eyevals(:,:,sacres.adjnid),3);
        for j = 1:size(pm,2)
            pm(:,j) = pm(:,j)-mean(pm(sid,j));
            nm(:,j) = nm(:,j)-mean(nm(sid,j));
        end
        hold on;
        plot(times,pm(:,usechan)-nm(:,usechan),'--');
        set(gca,'ylim',[-0.5 0.5]);
        if ~isempty(em.cptime)
            for j = 1:length(Expt.Trials)
                Expt.Trials(j).count = sum(Expt.Trials(j).Spikes > em.cptime(1) & Expt.Trials(j).Spikes < em.cptime(2));
            end
        end
        a = CalcCP([Expt.Trials(pid).count],[Expt.Trials(nid).count]);
        b = CalcCP([Expt.Trials(sacres.adjpid).count],[Expt.Trials(sacres.adjnid).count]);
        title(sprintf('Dashed lines = difference after removing %d/%d trials CP%.3f ->%.3f',length(sacres.removed),length(pid)+length(nid),a,b))
        end
    end
end

function whichsaccade(a,b,type, id, varargin)

    while ~isfigure(a)
        a = get(a,'Parent');
    end
    data = get(a,'UserData');
    
    tid = find([data.Trials.Trial] == data.trialid{type}(id));

    GetFigure('OneTrial');
    hold off;
    colors = 'rgbm';
    for j = 1:size(data.Trials(tid).EyeData,2)
    plot(data.times./10,data.Trials(tid).EyeData(:,j)-mean(data.Trials(tid).EyeData(1:100,j)),colors(j));
    hold on;
    end
    if isfield(data,'sacs')
    s = data.sacs{type}(id);
    t = s.start./10;
    fprintf('Sacc at %.3f in trial %d (R%.0f)\n',s.start./10000,data.trialid{type}(id), data.Trials(tid).RespDir);
    plot([t t],get(gca,'ylim'),'k:');
    end

function result =PlotSaccByChoice(Trials, Header, pidx, nidx, em, varargin)
   
result = [];
plotmode = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'mean',4)
        plotmode = 1;
    elseif strncmpi(varargin{j},'sumsac',4)
        plotmode = 2;
    elseif strncmpi(varargin{j},'cmpsac',4)
        plotmode = 4;
    elseif strncmpi(varargin{j},'cmppos',4)
        plotmode = 5;
    elseif strncmpi(varargin{j},'startpos',6)
        plotmode = 9;
    elseif strncmpi(varargin{j},'abspos',4)
        plotmode = 8;
    elseif strncmpi(varargin{j},'allverg',6)
        plotmode = 7;
    elseif strncmpi(varargin{j},'verg',4)
        plotmode = 6;
    end
    j = j+1;
end

sid = find(Header.emtimes > em.stattime(1)& Header.emtimes < em.stattime(2));
eid = find(Header.emtimes > em.stattime(3) & Header.emtimes < em.stattime(4));

npts = 50; %~70ms of data
for j = 1:length(Trials)
    T = Trials(j);
    for k = 1:length(T.Saccades)
        S = T.Saccades(k);
        [a,t] = min(abs(S.start-Header.emtimes));
        tid = max([1 t-npts]):t;
        hv(1) = mean(diff(T.EyeData(tid,1:2),[],2));
        vv(1) = mean(diff(T.EyeData(tid,3:4),[],2));
        [a,t] = min(abs(S.end-Header.emtimes));
        tid = t:min([t+npts length(Header.emtimes)]);
        vv(2) = mean(diff(T.EyeData(tid,3:4),[],2));
        hv(2) = mean(diff(T.EyeData(tid,1:2),[],2));
        Trials(j).Saccades(k).hv = diff(hv);
        Trials(j).Saccades(k).vv = diff(vv);
    end
end
if ismember(plotmode,[0 7]) %plot individual saccades
    psacs = [Trials(pidx).Saccades];
    id = find([psacs.start] < em.stattime(4) & [psacs.start] > em.stattime(1));
    [x,y] = pol2cart([psacs(id).dir], [psacs(id).size]);
    hold off;
    hw = range(y)/20;
    for j = 1:length(x)
        if plotmode == 7
            plot(x(j),psacs(id(j)).hv,'bo','buttondownfcn',{@whichsaccade, 1, id(j)});
        else
            plot(x(j),y(j),'o','buttondownfcn',{@whichsaccade, 1, id(j)});
        end
        hold on;
    end
    arrow([0 mean(x)],[0 mean(y)], 20, hw);
    n = 1;
    for j = pidx
        for k = 1:length(Trials(j).Saccades)
            data.trialid{1}(n) = Trials(j).Trial;
            n = n+1;
        end
    end


    nsacs = [Trials(nidx).Saccades];
    id = find([nsacs.start] < em.stattime(4) & [nsacs.start] > em.stattime(1));
    [x,y] = pol2cart([nsacs(id).dir], [nsacs(id).size]);
    for j = 1:length(x)
        if plotmode == 7
            plot(x(j),nsacs(id(j)).hv,'ro','buttondownfcn',{@whichsaccade, 2, id(j)});
        else
            plot(x(j),y(j),'ro','buttondownfcn',{@whichsaccade,2, id(j)});
        end
    end
    arrow([0 mean(x)],[0 mean(y)], 20, hw,'col','r');
    n = 1;
    for j = nidx
        for k = 1:length(Trials(j).Saccades)
            data.trialid{2}(n) = Trials(j).Trial;
            n = n+1;
        end
    end
    data.sacs{1} = psacs;
    data.sacs{2} = nsacs;
elseif ismember(plotmode,[2 4 5 6 6 7 8 9]) %plot sum of saccades for each trial
    hold off;
    for j = 1:length(Trials)
        hv = diff(Trials(j).EyeData(:,[1 2]),[],2);
        ch = mean(Trials(j).EyeData(:,em.hchan)',1);
        cv = mean(Trials(j).EyeData(:,em.vchan)',1);
        mdv(j) = mean(cv(eid)) - mean(cv(sid));
        mdh(j) = mean(ch(eid)) - mean(ch(sid));
        vpos(j,:) = [mean(cv(sid)) mean(cv(eid))];
        hpos(j,:) = [mean(ch(sid)) mean(ch(eid))];
        hverg(j,:) = [mean(hv(sid)) mean(hv(eid))];
    end

    for j = 1:length(pidx)
        S = Trials(pidx(j)).Saccades;
        if isempty(S)
            fprintf('Empty Saccades in Trial %d\n',pidx(j));
            x(j) = 0;
            y(j) = 0;
        else
        id = find([S.start] < em.stattime(4) & [S.start] > em.stattime(1));
        [a,b] = pol2cart([S(id).dir], [S(id).size]);
        x(j) = sum(a);
        y(j) = sum(b);
        hv(j) = sum([S(id).hv]);
        vv(j) = sum([S(id).vv]);
        if plotmode == 2
            plot(x(j),y(j),'o','buttondownfcn',{@whichsaccade, 1, j});
        elseif plotmode == 4
            plot(x(j),mdh(pidx(j)),'bo','buttondownfcn',{@whichsaccade, 1, j});
            plot(y(j),mdv(pidx(j)),'bs','buttondownfcn',{@whichsaccade, 1, j});
        elseif plotmode == 5
            plot(diff(hpos(pidx(j),:)),diff(vpos(pidx(j),:)),'bo','buttondownfcn',{@whichsaccade, 1, j});
        elseif plotmode == 8
            plot(hpos(pidx(j),1),hpos(pidx(j),2),'bo','buttondownfcn',{@whichsaccade, 1, j});
            plot(vpos(pidx(j),1),vpos(pidx(j),2),'bs','buttondownfcn',{@whichsaccade, 1, j});
        elseif plotmode == 9
            plot(hpos(pidx(j),1),vpos(pidx(j),1),'bo','buttondownfcn',{@whichsaccade, 1, j});
        elseif plotmode == 6
            plot(x(j),hv(j),'o','buttondownfcn',{@whichsaccade, 1, j});
        end
        hold on;
        end
    end
    allsac(pidx,:) = cat(1,x,y)';
    prefsac = [mean(x) mean(y)];
    hw = range(y)/20;
    arrow([0 mean(x)],[0 mean(y)], 20, hw);
    x = [];
    y = [];

    for j = 1:length(nidx)
        S = Trials(nidx(j)).Saccades;
        if isempty(S)
            fprintf('Empty Saccades in Trial %d\n',nidx(j));
            x(j) = 0;
            y(j) = 0;
        else
        id = find([S.start] < em.stattime(4) & [S.start] > em.stattime(1));
        [a,b] = pol2cart([S(id).dir], [S(id).size]);
        x(j) = sum(a);
        y(j) = sum(b);
        hv(j) = sum([S(id).hv]);
        vv(j) = sum([S(id).vv]);
        if plotmode == 4
        plot(x(j),mdh(nidx(j)),'ro','buttondownfcn',{@whichsaccade, 1, j});
        plot(y(j),mdv(nidx(j)),'rs','buttondownfcn',{@whichsaccade, 1, j});
        elseif plotmode == 8
            plot(hpos(nidx(j),1),hpos(nidx(j),2),'ro','buttondownfcn',{@whichsaccade, 2, j});
            plot(vpos(nidx(j),1),vpos(nidx(j),2),'rs','buttondownfcn',{@whichsaccade, 2, j});
        elseif plotmode == 5
            plot(diff(hpos(nidx(j),:)),diff(vpos(nidx(j),:)),'ro','buttondownfcn',{@whichsaccade, 2, j});
        elseif plotmode == 6
            plot(x(j),hv(j),'ro','buttondownfcn',{@whichsaccade, 2, j});
        elseif plotmode == 7
            plot(a,[S(id).hv],'ro','buttondownfcn',{@whichsaccade, 2, j});
        elseif plotmode == 9
            plot(hpos(nidx(j),1),vpos(nidx(j),1),'ro','buttondownfcn',{@whichsaccade, 1, j});
        else
        plot(x(j),y(j),'ro','buttondownfcn',{@whichsaccade, 2, j});
        end
        end
    end
    allsac(nidx,:) = cat(1,x,y)';
    if plotmode == 4
        allsac = cat(1,mdh,mdv)';
    elseif plotmode == 5 % 'cmppos' use pos diff to adjust
        allsac = cat(1,diff(hpos,[],2)',diff(vpos,[],2)')';
    end

    nullsac = mean(allsac(nidx,:));
    prefsac = mean(allsac(pidx,:));
    arrow([0 mean(x)],[0 mean(y)], 20, hw,'col','r','linewidth',2);
    diffsac = prefsac-nullsac;
    z = allsac * diffsac';
    crit = 100;
    diffmean = mean(z(pidx))-mean(z(nidx));
    while diffmean > 0
        b = prctile(z(pidx),100-crit);
        d = prctile(z(nidx),crit);
        nid = find(z(nidx) > b); % trials to include
        pid = find(z(pidx) < d);
        diffmean = mean(z(pidx(pid)))-mean(z(nidx(nid)));
        diffs(101-crit,1) = diffmean;
        crits(101-crit) = b;
        crit = crit-1;
    end
    if crit < 99
        crit = crit+1+diffs(end)./(diffs(end)-diffs(end-1));
    end
    b = prctile(z(pidx),100-crit);
    d = prctile(z(nidx),crit);
    nid = find(z(nidx) > b); % trials to include
    pid = find(z(pidx) < d);
    result.adjpid = pidx(pid);
    result.adjnid = nidx(nid);
    result.removed = [setdiff(pidx, pidx(pid)) setdiff(nidx, nidx(nid))];
    result.diffmean(1) = mean(z(pidx(pid)))-mean(z(nidx(nid)));
    result.diffmean(2) = mean(diff(hpos(pidx(pid),:),[],2))-mean(diff(hpos(nidx(nid),:),[],2));
    result.diffmean(3) = mean(diff(vpos(pidx(pid),:),[],2))-mean(diff(vpos(nidx(nid),:),[],2));
    result.diffmean(4) = mean(diff(hverg(pidx(pid),:),[],2))-mean(diff(hverg(nidx(nid),:),[],2));
    result.ids = [Trials.id];
    
    pm = (mean(cat(3,Trials(result.adjpid).EyeData),3));
    nm = (mean(cat(3,Trials(result.adjnid).EyeData),3));

    if plotmode == 2
        for j = 1:length(nid)
            plot(allsac(nidx(nid(j)),1),allsac(nidx(nid(j)),2),'ro',...
                'markerfacecolor','r','buttondownfcn',{@whichsaccade, 2, j});
        end
        for j = 1:length(pid)
            plot(allsac(pidx(pid(j)),1),allsac(pidx(pid(j)),2),'bo',...
                'markerfacecolor','b','buttondownfcn',{@whichsaccade, 2, j});
        end
        title('Removing Open trials eliminates differences in sum(saccades)');
    elseif plotmode == 6
        ylabel('summed Saccadic Change in vergence');
        xlabel('summed Saccadic Change in conjug');
    elseif plotmode == 7
        ylabel('Saccadic Change in vergence');
        xlabel('Saccadic Change in conjug');
        title('Each circle is one saccade');
    end
    data.trialid{1} = [Trials(pidx).Trial];
    data.trialid{2} = [Trials(nidx).Trial];
elseif plotmode == 1 %calculate change in position over trial
    sid = find(Header.emtimes > em.stattime(1) & Header.emtimes < eb.stattime(3));
    eid = find(Header.emtimes > em.stattime(3) & Header.emtimes < em.stattime(4));
    aid = find(Header.emtimes > em.stattime(1) & Header.emtimes < em.stattime(4));
    hold off;
    starts = [];
    sizes = [];
    for k = 1:length(pidx)
        j = pidx(k);
        if isempty(Trials(j).Saccades)
            aid = [];
        else
            starts = [starts [Trials(j).Saccades.start]];
            sizes = [sizes [Trials(j).Saccades.size]];
            aid = find([Trials(j).Saccades.start] > em.stattime(1) & [Trials(j).Saccades.start] < em.stattime(4));
        end
        if length(aid)
        [x,y] = pol2cart([Trials(j).Saccades(aid).dir],[Trials(j).Saccades(aid).size]);
        sdv(j) = sum(y);
        else
            sdv = NaN;
        end
        ch = mean(Trials(j).EyeData(:,em.hchan)',1);
        cv = mean(Trials(j).EyeData(:,em.vchan)',1);
        mdv(j) = mean(cv(eid)) - mean(cv(sid));
        mdh(j) = mean(ch(eid)) - mean(ch(sid));
        if plotmode == 3 %compare mean end -start vs sum of saccades
            plot(mdv(j),sdv(j),'o','buttondownfcn',{@whichsaccade,1, k});
        else
            plot(mdh(j),mdv(j),'o','buttondownfcn',{@whichsaccade,1, k});
        end
        hold on;
    end
    for k = 1:length(nidx)
        j = nidx(k);
        if isempty(Trials(j).Saccades)
            aid = [];
        else
            aid = find([Trials(j).Saccades.start] > em.stattime(1) & [Trials(j).Saccades.start] < em.stattime(4));
        end
        if length(aid);
            [x,y] = pol2cart([Trials(j).Saccades(aid).dir],[Trials(j).Saccades(aid).size]);

            sdv(j) = sum(y);
        else
            sdv(j) = NaN;
        end
        ch = mean(Trials(j).EyeData(:,em.hchan)',1);
        cv = mean(Trials(j).EyeData(:,em.vchan)',1);
        mdv(j) = mean(cv(eid)) - mean(cv(sid));
        mdh(j) = mean(ch(eid)) - mean(ch(sid));
        if plotmode == 3 %compare mean end -start vs sum of saccades
            plot(mdv(j),sdv(j),'ro','buttondownfcn',{@whichsaccade,2, k});
        else
            plot(mdh(j),mdv(j),'ro','buttondownfcn',{@whichsaccade,2, k});
        end
        hold on;
    end
    plot(mean(mdh(pidx)),mean(mdv(pidx)),'s','markerfacecolor','b');
    plot(mean(mdh(nidx)),mean(mdv(nidx)),'s','markerfacecolor','r');
    sdir = [mean(mdh(pidx))-mean(mdh(nidx)) mean(mdv(pidx))-mean(mdv(nidx))]; 
    shifts = cat(2, mdh, mdv);
    data.trialid{1} = [Trials(pidx).Trial];
    data.trialid{2} = [Trials(nidx).Trial];
        
end
 data.Trials = Trials;
 data.times = Header.emtimes;
 set(gcf,'UserData',data);
 axis('image');
 

 
 
    function Time2Sample(t, TimeList)
        
        
        
        