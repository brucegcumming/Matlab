function PlotAllExpts(Expts, xc)
%PlotAllExps plots results that combine multiple expts. 
%Currently only plots cross correlograms

PlotXcorrs(Expts, xc, 'xcorr-test');

function PlotXcorrs(Expts, xc, plottype)

    hold off;
    for j = 1:length(xc)
        Ea = Expts{xc(j).cells(1)};
        Eb = Expts{xc(j).cells(2)};
        xcs(j) = xc(j).rsc(1);
        if strcmp(plottype,'xcorr-distance')
        plot(abs(xc(j).probesep),xc(j).rsc(1),'o','buttondownfcn',{@HitXcorr, j});
        elseif strcmp(plottype,'xcorr-cp')
        plot((0.5 - Eb.cp)*(0.5-Ea.cp),xc(j).rsc(1),'o','buttondownfcn',{@HitXcorr, j});
        elseif strcmp(plottype,'xcorr-test')
        plot(abs(Eb.cp-Ea.cp),xc(j).rsc(1),'o','buttondownfcn',{@HitXcorr, j});
        end
        hold on;
    end
        PredictChoices(Expts,xc);
    title(sprintf('Mean %.3f',mean(xcs(~isnan(xcs)))));

   
   function PredictChoices(Expts, xc);
    ids = [];
    for j = 1:length(Expts)
        zid = find([Expts{j}.Trials.ob] > 120 & abs([Expts{j}.Trials.RespDir]) == 1);
        zids{j} = [Expts{j}.Trials(zid).id];
        ids = cat(2,ids,zids{j});
        cps(j) = Expts{j}.cp;
    end
    [a,b] = Counts(ids);
    bar(b,a);
    c = max(a);
    smw = 10;
    for j = 1:floor(c/2);
    smw = 10;
    t = smooth(a >= c,smw);
    while max(t) > 0.999
        smw = smw+1;
        t = conv(a >= (c-j),ones(1,smw)./smw);
    end
    maxlen(j) = smw-1;
    end
    id = find(maxlen > max(maxlen)/2);
    nc = id(1);
    nt = maxlen(id(1));
    t = conv(a >= (c-j),ones(1,smw));
    w = ceil(smw/2);
    [a,c] = max(t(w:end));
    w = floor(nt/2)-1;
    blist = b([c-w:c+w]);
    allids = [];
    for j = 1:length(zids)
        if sum(ismember(zids{j},blist)) > nt*0.8 
            includecell(j) = 1;
            if isempty(allids)
                allids = zids{j};
            else
                allids = intersect(zids{j},allids);
            end
        end
    end
length(allids);    
includecell = find(includecell);            
for j = 1:length(includecell)
    c = includecell(j);
    id = find(ismember([Expts{c}.Trials.id],allids));
    counts(:,j) = [Expts{c}.Trials(id).count];
    choices(:,j) = [Expts{c}.Trials(id).RespDir];
end

wts = cps(includecell) -0.5;
cp = WeightedCP(wts, counts,choices(:,1));
fwts = fminsearch( @WeightedCP, wts, optimset , counts, choices(:,1))
x = [-1 -0.5 -0.25 0 0.5 0.5 1];
[a,b,c,d,e] = ndgrid(x,x,x,x,x);
for j = 1:length(a(:))
fwts = fminsearch( @WeightedCP, [a(j) b(j) c(j) d(j) e(j)], optimset , counts, choices(:,1));
cpi(j) = WeightedCP(fwts, counts,choices(:,1));
end
cp = WeightedCP(fwts, counts,choices(:,1))
[cp,j] = min(cpi);
fwts = fminsearch( @WeightedCP, [a(j) b(j) c(j) d(j) e(j)], optimset , counts, choices(:,1))
cp = WeightedCP(fwts, counts,choices(:,1))
[sortcp, id] = sort(abs(cps(includecell)-0.5),'descend');
for j = 4:length(id)
t = treefit(counts(:,id(1:j)),choices(:,1),'method','classification');
x = eval(t,counts(:,id(1:j)));
score = length(strmatch('1',x(choices(:,1)==1))) + length(strmatch('-1',x(choices(:,1)==-1)));
score = score./size(choices,1)
   end
[sortcp, id] = sort(abs(cps(includecell)-0.5),'descend');
for j = 2:length(id)
    C = classify(counts(:,id(1:j)),counts(:,id(1:j)),choices(:,1),'quadratic');
    scores = sum(C == choices(:,j))./size(choices,1)
end
C = classify(counts,counts,choices(:,1),'quadratic');

function cp = WeightedCP(w, counts,choices)
    
    sumcount = counts * w';
    cp = CalcCP(sumcount(choices == 1),sumcount(choices == -1));
        
    

