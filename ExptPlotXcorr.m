function ExptPlotXcorr(xcorrs,bycell)

callback = @HitXCorrAll;

    ClearPlot;
    X.crit.synci = 1;
    X.plot.synci = 0;
    X.plot.xcorrlabeltop = 0;
    X.corrs = xcorrs;
    if bycell
%        xcorrs = xcorrs([xcorrs.valid] == 1);
        cellids = cat(1,xcorrs.cells);
        cellids(isnan(cellids)) = 0;
    else
        cellids = cat(1,xcorrs.probes);
    end
%    exids = cat(1,xcorrs.eid);
    exids = ones(size(cellids,1),1);
    expts = 1;
    weights = prod(cat(1,xcorrs.ntrials)');
    weights = cat(1,xcorrs.ntrials);
    cells = unique(cellids);
    cells = cells(cells > 0);
    probes = cat(1,xcorrs.probes);
    np = length(cells);
    ns = 0;
    for j = 1:length(cells)
        ida = find(cellids(:,1) == cells(j) & ismember(exids,expts));
        idb = find(cellids(:,2) == cells(j) & ismember(exids,expts));
        cellpos(j) = (sum(probes(ida,1))+sum(probes(idb,2)))./(length(ida)+length(idb));
        for k = 1:j
            ida = find(cellids(:,1) == cells(j) & cellids(:,2) == cells(k) & ismember(exids,expts));
            idb = find(cellids(:,2) == cells(j) & cellids(:,1) == cells(k) & ismember(exids,expts));
            if length(ida)+length(idb) > 0
                separation(j,k) =( sum([xcorrs(ida).probesep]) - sum([xcorrs(idb).probesep]))./(length(ida)+length(idb));
            else
                separation(j,k) = cellpos(j)-cellpos(k);
            end
            separation(k,j) = -separation(j,k);
        end
    end
    order = sum(separation < 0);
    [a,b] = sort(cellpos);
    icells = cells; %unsorted
    cells = cells(b);
    cellpos = cellpos(b);
    ta = ceil(length(xcorrs(1).xc)/2);
    tb = ta + length(xcorrs(1).xc)-1;
    midpt = 1+floor(length(xcorrs(1).xc)/2);
    for j = 1:length(cells)
        for k = 1:j
            id = find((cellids(:,1) == cells(j) & cellids(:,2) == cells(k)) | (cellids(:,2) == cells(j) & cellids(:,1) == cells(k)));
            id = id(ismember(exids(id),expts));
            if length(id)
                if length(id) > 1
                    xc = WeightedSum(cat(1,xcorrs(id).xc),weights(id));
                else
                    xc = xcorrs(id).xc;
                end
                if j == k
                    xc(midpt) = 0;
                end
                mysubplot(np,np,k+(j-1)*np,'leftmargin',0.02);
                h = plot(ta:tb,xc,'k-','linewidth',2);
                synci = SyncIndices(xc);
                if synci(2) < X.crit.synci
                    set(h,'color','r');
                end
                if ~isnan(synci(2))
                    ns = ns+1;
                    syncis(ns,1:length(synci)) = synci;
                    syncis(ns,3) = cells(j);
                    syncis(ns,4) = cells(k);
                end
                axis('tight');
                xl = get(gca,'xlim');
                yl = get(gca,'ylim');
                set(gca,'xtick',[],'ytick',[],'buttondownfcn',{callback,bycell,[cells(j) cells(k)], expts});
                set(h,'buttondownfcn',{callback,bycell,[cells(j) cells(k)],expts});
                if k == j
                    if j == 1 || X.plot.xcorrlabeltop == 1
                        title(sprintf('Cell%d',cells(k)));
                    else
                        ii = find(icells == cells(k-1));
                        ij  = find(icells == cells(j));
                        if bycell
                            h = text(xl(1),yl(2),sprintf('C%d at %.1f',cells(k),cellpos(k)));
                            set(h,'HorizontalAlignment','left','verticalalignment','bottom');
                        else
                            title(sprintf('P%d',cells(k)));
                        end
                    end
                end
                if k == 1
                    ylabel(sprintf('Cell%d',cells(j)));
                end
            end
        end
    end
    if X.plot.synci
    mysubplot(2,2,2);
    myscatter(syncis(:,1),syncis(:,2),'o','ids',syncis(:,3:4));
    end
    X.toplevel = gcf;
    set(X.toplevel,'UserData',X);

function synci = SyncIndices(xc)
midpt = ceil(length(xc)./2);
synci(1) = mean([xc(midpt-1) xc(midpt+1)])./xc(midpt);
synci(2) = mean([xc(1:midpt-10) xc(midpt+10:end)])./xc(midpt);


function HitXCorrAll(a,b, bycell, cells, expts)
X = GetDataFromFig(a);
mysubplot(2,2,2);
cellids = cat(1,X.corrs.cells);
aid = find(cellids(:,1) == cells(1) & cellids(:,2) == cells(2));
bid = find(cellids(:,1) == cells(2) & cellids(:,2) == cells(1));
if isempty(bid)
    id = aid;
else
    id = bid;
end
if cells(1) == cells(2)
    midpt = 1+floor(length(X.corrs(id).xc)/2);
    X.corrs(id).xc(midpt) = 0;
end
plot(X.corrs(id).xc);
axis('tight');
set(gca,'xtick',[],'ytick',[]);
