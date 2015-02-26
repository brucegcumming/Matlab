function CompareAutoClusters(DATA, expts)


%X.efficacy is a NclManuaul x 2 * NclAuto array
%X.matches(k,2) is the cluser in Cb that matches Cluster k in Ca
%X.matches(k,1) is the efficacy
matches = [];
for e = 1:length(expts)
    ei = expts(e);
    ci = find(ismember(DATA.exptid,ei));
    Ca = PC.CheckClusterLoaded(DATA,ci,'auto');
    C = PC.CheckClusterLoaded(DATA,ci);
    for j = 1:length(C{e});
        X{e,j} = PC.CompareClusters(C{ci}{j},Ca{ci}{j},DATA, 'noplot');
        X{e,j}.exptno = ei;
        X{e,j}.probe = Ca{e}{j}.probe(1);
        dp(ci,j) = PC.DistanceMeasure(C{ci}{j},1,DATA.mahaltype);
        if isempty(X{e,j}.efficacy)
            errors(ci,j) = 1;
            mycl = 1;
        else
            matches(j,1:size(X{e,j}.match,1)) = X{e,j}.match(:,1);
            matchid(j,1:size(X{e,j}.match,1)) = X{e,j}.match(:,2);
            mycl = X{e,j}.match(1,2); %cl num of manual cluster that matches #1
        end
        dp(ci,j) = PC.DistanceMeasure(C{ci}{j},mycl,DATA.mahaltype);
        adp(e,j) = PC.DistanceMeasure(Ca{ci}{j},mycl,DATA.mahaltype);
        automap(e,j) = C{e}{j}.auto;
    end
    
    np = size(matches,1);
    for j = 1:size(matches,2)
        allmatches(e,1:np,j) = matches(:,j);
        allmatchid(e,1:np,j) = matchid(:,j);
    end
end

setappdata(DATA.toplevel,'automatches',allmatches);
GetFigure(DATA.tag.all,'parent',DATA.toplevel);
imagesc(squeeze(allmatches(:,:,1)));
F = GetFigure(DATA.tag.popscatter,'parent',DATA.toplevel);
hm = findobj(F,'type','uimenu','tag','ComparePlotMenu');
if isempty(hm)
    hm = uimenu(F,'Label','Plot','tag','ComparePlotMenu');
    uimenu(hm,'Label', 'Efficacies','callback',{@ReplotComparison, 'efficacies'});
    uimenu(hm,'Label', 'Cell Match Efficacies','callback',{@ReplotComparison, 'cellefficacies'});
    uimenu(hm,'Label', 'Cell All Efficacies','callback',{@ReplotComparison, 'allcellefficacies'});
    uimenu(hm,'Label', 'Matches','callback',{@ReplotComparison, 'matches'});
    uimenu(hm,'Label', 'Isolation','callback',{@ReplotComparison, 'isolation'});
end


PlotData.X = dp;
PlotData.Y = allmatches(:,:,1);
PlotData.allmatches = allmatches;
PlotData.allmatchid = allmatchid;
PlotData.CompData = X;

setappdata(F,'PlotData',PlotData);
plot(dp,allmatches(:,:,1),'o','buttondownfcn',{@HitScatter});

function ReplotComparison(a,b, type)

F = GetFigure(a);
DATA = GetDataFromFig(F);

PlotData = getappdata(gcf,'PlotData');
PlotData.plottype = type;
X = PlotData.CompData;
typelabels = {'Unmatched Cell' 'Cell' 'Not cell' 'Nearest' 'Best' 'Best(MU'};

if sum(strcmp(PlotData.plottype,{'efficacies' 'cellefficacies'}))
    allefficacy = [];
    for j = 1:size(X,1)
        e = find(ismember(DATA.exptid,X{j}.exptno));
        ce = find(ismember(DATA.CellDetails.exptids,X{j}.exptno));
        p = X{j}.probe;
        for a = 1:size(X{j}.efficacy,1)
            [abest, abestid]  = max(squeeze(X{j}.efficacy(a,2,:)));
            [best, bestid]  = max(squeeze(X{j}.efficacy(a,1,:)));
            if PC.isacell(DATA, ce,p,a)
                if best < 0.8 && abest < 0.8
                    result = 1;
                else
                    result = 2;
                end
            else
                result = 3;
            end
            for b = 1:size(X{j}.efficacy,3)
                allefficacy(end+1,1:2) = X{j}.efficacy(a,:,b);
                allefficacy(end,3) = a;
                allefficacy(end,4) = b;
                allefficacy(end,5) = X{j}.probe;
                allefficacy(end,6) = X{j}.exptno;
                matches(size(allefficacy,1)) = result;
                if b == bestid
                    matches(size(allefficacy,1)) = result+3;
                end
            end
        end
    end
    colors = mycolors('white');
    PlotData.allefficacy = allefficacy;
    PlotData.X = allefficacy(:,1);
    PlotData.Y = allefficacy(:,2);
    types = unique(matches);
    if strcmp(PlotData.plottype,'cellefficacies')
        types = [4 5];
    end
    hold off;
    for j = 1:length(types)
        id = find(matches == types(j));
        labels{j} = typelabels{types(j)};
        plot(allefficacy(id,1),allefficacy(id,2),'o','buttondownfcn',@HitScatter,'color',colors{types(j)});
        hold on;
    end
    legend(labels);
    setappdata(gcf,'PlotData',PlotData);
elseif strcmp(type, 'allcellefficacies')
    nc = 0;
    nclusters = ones(1,DATA.nprobes);
    c = squeeze(sum(abs(DATA.CellList),1));
    c(:,1) =1;
    cols = cumsum([1 sum(c' > 0)]);
    matchim = zeros(size(X,1),max(cols));
    for ei = 1:size(X,1)
    for j = 1:size(X,2)
        A = X{ei,j};
            e = find(ismember(DATA.exptid,A.exptno));
            ce = find(ismember(DATA.CellDetails.exptids,A.exptno));
            p = A.probe;
            [a,b,clid] = PC.isacell(DATA,ce,p);
            clid = clid(clid > 0);
            if max(clid) > nclusters(p)
                nclusters(p) = max(clid);
            end
            c = cols(p);
            for k = 1:length(clid)
                nc = nc+1;
                c = cols(p) + clid(k)-1;
                id = find(A.mlst == clid(k)+1);
                if size(A.efficacy,1) < id
                    eff = [];
                else
                    eff = squeeze(A.efficacy(id,:,:));
                end
                if size(eff,1) < 2
                    if isempty(eff)
                        Im(e,c) = 3;
                    else
                        Im(e,c) = 2;
                    end
                else
%most of the auto must be in the manual - if the auto is a much bigger
%set of events, its not a match.  The other way round is less clear
%an auto cut may subdividea a manual SD
                    [a(1),b] = max(eff(2,:)); %Fraction of auto cut htat are cell
                    [a(2),d] = max(eff(1,:));  %fraction of cell in autocut
                    allefficacy(nc,1) = eff(1,b);
                    allefficacy(nc,2) = eff(2,b);
                    if (a(1) > 0.8 && a(2) > 0.2) %if v small # of cell is is auto,
                        Im(e,c) = 1;
                    else
                        Im(e,c) = 2;
                    end
                    matchim(e,c) = b;
                end
            end
            pid = cols(p):cols(p+1)-1;
            [n, cls] = Counts(matchim(e,pid));
            if max(n) > 1 %same auto matches two cells
                id = find(n > 1 & cls > 0);
                for k = 1:length(id)
                    cid = find(matchim(e,pid) == cls(id(k)));
                    Im(e,pid(cid)) = 4;
                end
            end
    end
    end
    PlotData.X = allefficacy(:,1);
    PlotData.Y = allefficacy(:,2);
    PlotData.matchim = matchim;
    PlotData.cols = cols;
    hold off;
    plot(allefficacy(:,1),allefficacy(:,2),'o','buttondownfcn',@HitScatter);
    GetFigure(DATA.tag.all,'parent',DATA.toplevel);
    hold off;
    h = imagesc(Im);
    for j = 2:length(cols)
        line([cols(j)-0.5 cols(j)-0.5],[0 size(Im,2)],'color','w');
    end
            
            
    set(h,'buttondownfcn',{@HitImage,'cellmatches'});
    setappdata(gcf,'PlotData',PlotData);
end

function HitImage(a,b, type)

PlotData = getappdata(gcf,'PlotData');
DATA = GetDataFromFig(gcf);
xy = get(gca,'currentpoint');
ex = round(xy(1,2));
col = round(xy(1,1));
id = find(PlotData.cols <= col);
p = PlotData.cols(id(end));
probe = id(end);
clid(1) = 1+col-p;
clid(2) = PlotData.matchim(ex,col);
PlotAutoCompare(DATA,PlotData,ex,probe,clid);

function PlotAutoCompare(DATA,PlotData,e,p, cls)

AutoClusters = getappdata(DATA.toplevel,'AutoClusters');
autoC = AutoClusters{e}{p};
Clusters = getappdata(DATA.toplevel,'Clusters');
C = Clusters{e}{p};
DATA.plot.autoxy = 0;
PC.SetFigure(DATA,DATA.tag.xyplot,'front');
DATA.currentpoint = [e p];
DATA = PC.ConditionalPlotXY(DATA, C, 0,'cluster',cls(1));
PC.SetFigure(DATA,DATA.tag.autoxyplot,'front');
hold off;

PC.PlotClusterXY(DATA,autoC,'cellid',cls(2)-1);
A = PlotData.CompData{e,p};
cellno = DATA.CellList(e,p,cls(1));
id = find(A.mlst == cls(1)+1);
if cls(2) == 0
    [a,b] = max(A.efficacy);
    if size(A.efficacy,3) == 1
        fprintf('E%dP%dCl%d Cell%d AutoCut Failed = all one codeEff %.2f %.2f\n',e,p,cls(1),cellno,A.efficacy(id,1,1),A.efficacy(id,2,1));
    else
        if size(A.efficacy,1) < id
            fprintf('Nonmatch E%dP%dCl%d Cell%d Missing efficacy\n',e,p,cls(1),cellno);
            cid = unique(AutoClusters{e}{p}.clst);
        else
        [a,b] = max(max(squeeze(A.efficacy(id(1),:,:))));
        fprintf('Nonmatch E%dP%dCl%d Cell%d->AutoBest%d %.2f A->C %.2f\n',e,p,cls(1),cellno,b(1),A.efficacy(id,1,b(1)),A.efficacy(id,2,b(1)));
        end
    end
else
    nc = [sum(C.clst == 1+cls(1)) sum(autoC.clst ==cls(2))];
    fprintf('E%dP%dCl%d Cell%d(%d)->Auto%d(%d) %.2f A->C %.2f\n',e,p,cls(1),cellno,nc(1),cls(2),nc(2),A.efficacy(id,1,cls(2)),A.efficacy(id,2,cls(2)));
end

function HitScatter(a,b)

F = GetFigure(a);
DATA = GetDataFromFig(F);
pos = get(gca,'currentpoint');
xy = pos(1,1) + i* pos(1,2);
P = getappdata(F,'PlotData');
D = abs(P.X+i*P.Y - xy);
[d, id] = min(D(:));
if sum(strcmp(P.plottype,{'efficacies' 'cellefficacies'}))
    x = P.allefficacy(id,:);
    expt = x(6);
    p = x(5);
    fprintf('E%dP%d Cluster %d, Auto %d\n',expt,p,x(3),x(4));
    e = find(ismember(DATA.exptid,expt));
else
    [e,p] = ind2sub([size(P.allmatches,1) size(P.allmatches,2)],id);
    fprintf('E%dP%d\n',e,p);
    x(3) = P.allmatches(id,1);
end
DATA = GetDataFromFig(F);
AutoClusters = getappdata(DATA.toplevel,'AutoClusters');
autoC = AutoClusters{e}{p};
Clusters = getappdata(DATA.toplevel,'Clusters');
C = Clusters{e}{p};
DATA.plot.autoxy = 0;
DATA = PC.ConditionalPlotXY(DATA, C, 0,'cluster',x(3));
PC.SetFigure(DATA,DATA.tag.autoxyplot,'front');
hold off;
PC.PlotClusterXY(DATA,autoC,'cellid',x(4)-1);

%DATA = PC.ConditionalPlotXY(DATA, autoC, 0,'cluster',x(4)-1);

%PC.ShowData(DATA, e,p);
F = GetFigure('CompareAuto');
PC.CompareClusters(C,autoC, DATA,'clusters', [x(3) x(4)-1]);
CompareAutoMeans(DATA,e,p,x(3),x(4));

function CompareAutoMeans(DATA, e, p, ca, cb)

Clusters = getappdata(DATA.toplevel,'Clusters');
AutoClusters = getappdata(DATA.toplevel,'AutoClusters');
Ca = PC.GetClusterInfo(Clusters{e}{p},ca);
Cb = PC.GetClusterInfo(AutoClusters{e}{p},cb-1);
GetFigure(DATA.tag.spkmean);
mysubplot(2,1,1);
plot(Ca.MeanSpike.ms');
title(sprintf('Cluster %d',ca),'verticalalignment','top');
mysubplot(2,1,2);
plot(Cb.MeanSpike.ms');
title(sprintf('Auto Cluster %d',cb),'verticalalignment','top');

