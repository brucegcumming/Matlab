function X = CompareClusters(Ca, Cb, DATA, varargin)
%X = CompareClusters(Ca, Cb, DATA.. Compare two Clusters in PlotClusters.
% (e.g. auto and manual cut
%Calculates Efficacy to check for identity match
%X.efficacy is a Ncla x 2 * Nclb array
%X.matches(k,2) is the cluser in Cb that matches Cluster k in Ca
%X.matches(k,1) is the efficacy
%always finds match - just the cluster with the highest efficacy
plottype = 'all';
includehash = 2;
useclusters = [];
X = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'cluster',5)
        j = j+1;
        useclusters = varargin{j};
    elseif strncmpi(varargin{j},'noplot',5)
        plottype = 'none';
    end
    j = j+1;
end

if isfield(Ca,'clst') && isfield(Cb,'clst')
    if isfield(Ca,'t') &&  length(Ca.t) == length(Ca.clst)
        Ca.times = Ca.t;
    end
    if isfield(Cb,'t') &&  length(Cb.t) == length(Cb.clst)
        Cb.times = Cb.t;
    end
    if length(Ca.times) == length(Ca.clst) %ClusterDetails
        if isfield(Ca,'exptid')
            e = Ca.exptid;
        elseif isfield(Cb,'exptid')
            e = Cb.exptid;
        else
            e = find(ismember(DATA.exptid,Ca.exptno));
        end
        p = Ca.probe(1);
        nac = unique(Cb.clst);
        Cborder = nac;
        nc = unique(Ca.clst);
        Caorder = nc;
        if includehash == 2 %only fro Cb
            nc = nc(nc > 1);
        elseif includehash == 0
            nac = nac(nac>1);
            nc = nc(nc > 1);
        end
        if ~isempty(useclusters)
            nc = 1+useclusters(1);
            Caorder(nc) = Caorder(end);
            Caorder(end) = nc;
        end
        if length(useclusters) > 1
            nac = 1+useclusters(2);
            Cborder(nac) = Cborder(end);
            Cborder(end) = nac;
        end
        if isempty(nac) || isempty(nc)
            X.efficacy = [];
            X.match = [];
        else
            for k = 1:length(nac)
                tb = Cb.times(Cb.clst ==nac(k));
                for j = 1:length(nc)
                    ta = Ca.times(Ca.clst ==nc(j));
                    [xc(j,:), D] = xcorrtimes(ta, tb);
                    efficacy(j,:,k) = D.efficacy;
                end
                [a,b] = max(max(efficacy(:,:,k),[],2));
                fprintf('E%dP%d Cl%d->%d Eff%.3f,%.3f\n',Ca.exptno,Ca.probe(1),nac(k)-1,nc(b)-1,a,min(efficacy(b,:,k)));
                X.match(k,:) = [a b];
            end
            X.efficacy = efficacy;
            X.alst = nac; %list of auto clusters
            X.mlst = nc; %list of manual clusters
        end
        if strcmp(plottype,'all')
            subplot(2,2,1);
            plot(xc');
            [a,b] = max(max(efficacy(:,:,1)'));
            title(sprintf('Eff%.3f,%.3f',a,min(efficacy(b,:,1))));
            subplot(2,2,2);
            AllSpikes = PC.CheckAllSpikes(DATA, e, p, 'allprobes');
            PC.QuickSpikes(DATA, AllSpikes{e,p}, Ca,'oneprobe','order',Caorder);
            nev = length(Ca.clst);
            if ~isempty(useclusters)
                title(sprintf('Cluster%d %d/%d spikes',useclusters(1),sum(Ca.clst == useclusters(1)+1),nev));
            else
                title(sprintf('Spikes %s/%d',sprintf(' %d',Counts(Ca.clst)),nev));
            end
            aname = strrep(AllSpikes{e,p}.filename,'/Spikes/','/AutoSpikes/');
            Spikes = ReadSpikeFile(aname);
            subplot(2,2,3);
            PC.QuickSpikes(DATA, Spikes, Cb,'order', Cborder);
            nev = length(Cb.clst);
            if length(useclusters) > 1
                title(sprintf('AutoCut Cl%d %d/%d spikes',useclusters(2),sum(Cb.clst == useclusters(2)+1),nev));
            else
                title(sprintf('AutoCut %d/%s',length(Spikes.values),nev));
            end
            subplot(2,2,4);
            hist(difftimes(ta,tb),[0:0.02:0.5]);
            axis tight;
            GetFigure('CompareTimes','parent',DATA.toplevel)
            nc = unique(Cb.clst);
            nc = nc(nc > 0);
            hold off;
            for j = 1:length(nc)
                id = find(Cb.clst == nc(j));
                plot(Cb.times(id),Cb.triggerV(id),'.','color',DATA.colors{nc(j)});
                hold on;
            end
            if isfield(Ca,'triggerV')
            nc = unique(Ca.clst);
            nc = nc(nc > 0);
            for j = 1:length(nc)
                id = find(Ca.clst == nc(j));
                plot(Ca.times(id),Ca.triggerV(id),'o','color',DATA.colors{nc(j)});
            end
            end
            set(gca,'buttondownfcn',{@ReportTrial});
        end
    elseif isfield(Ca,'times') && isfield(Cb,'times')
        [a,b] = CalcEfficacy(Ca.times, Cb.times);
        X.efficacy = a;
    else
        X.efficacy = [];
    end
else
    X.efficacy = [];
end

function ReportTrial(a,b)
DATA = GetDataFromFig(a);
x = get(gca,'currentpoint');
t = x(1,1);
for j = 1:length(DATA.Expt.Trials)
    Starts(j) = DATA.Expt.Trials(j).Start(1);
    Ends(j) = DATA.Expt.Trials(j).End(end);
end
    
sid = find(Starts < t .*10000);
eid = find(Ends > t .*10000);

if ~isempty(sid) && ~isempty(eid)
    tid = DATA.Expt.Trials(sid(end)).id;
    Tn = DATA.Expt.Trials(sid(end)).Trial;
    fprintf('Trials %d Id%d Trial %d at %.2f\n',sid(end),tid,Tn,t)
end


function diffs = difftimes(ta, tb)

for j = 1:length(ta)
    diffs(j) = min(abs(ta(j)-tb));
end
