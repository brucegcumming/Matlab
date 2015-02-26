function [cells, xcorrs] = PlotAllXCorr(DATA, DataClusters, cells, varargin)
callback = @PlotXcorr;
j = 1;
ids = [];
while j <= length(varargin)
    if strcmp(varargin{j},'callback')
        j = j+1;
        callback = varargin{j};
    elseif strcmp(varargin{j},'sublist')
        j = j+1;
        ids = varargin{j};
    end
    j = j+1;
end
if isempty(ids)
    ids = 1:length(cells);
end
    ClearPlot;
    setappdata(DATA.toplevel,'xcCellList',cells);
%    xpts = linspace(DATA.spts(1),DATA.spts(end),401);
    nc = length(cells);
    nxc = 0;
    for j = 1:length(cells)
    for k = 1:j
        if ismember(k,ids) || ismember(j,ids)
        mysubplot(nc,nc,(j-1) * nc+k);
%        xoff = floor((j-1)/6) .* length(DATA.spts);
%        yoff = rem(j-1,6) .* DATA.vsep;

        P = DataClusters{cells(j).p};
        if cells(j).cl > 1
            P  = P.next{cells(j).cl-1};
        end
        Q = DataClusters{cells(k).p};
        if cells(k).cl > 1
            Q  = Q.next{cells(k).cl-1};
        end
        [xc, details] = xcorrtimes(P.times,Q.times);
        xpts = details.xpts;
        nxc = nxc+1;
        xcorrs(nxc).shapexc = xShapeCorr(P,Q);
        xcorrs(nxc).efficacy = details.efficacy;
        if (j == k)
            xc(details.midpt) = 0;
        end
        xcorrs(nxc).xc = xc;
        xcorrs(nxc).p = [j k];
        h = plot(xpts,xc,'k','linewidth',2);
        [a,b]= max(xc);
        if max(details.efficacy) > DATA.crit.synci(2)
            set(h,'color','r');
        end
        set(gca,'xtick',[],'ytick',[],'buttondownfcn',{callback, j,k});
        set(gca,'UserData',[j k]);
        set(h,'buttondownfcn',{callback, j,k});
        axis('tight');
        drawnow;
        if k == j 
            if isfield(P,'fitdprime') && P.fitdprime(1) < -2
                set(h,'color','r');
            else
                set(h,'color','k');
            end
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            h = text(xl(2),yl(2),sprintf('%d/%d',cells(j).p,cells(j).cl),...
                'horizontalalignment','right',...
                'verticalalignment','bottom');
            if j == 1
                set(h,'horizontalalignment','right',...
                'verticalalignment','bottom');
            end
        elseif k == max(ids)
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            h = text(xl(2),yl(2),sprintf('%d/%d',cells(j).p,cells(j).cl),...
                'horizontalalignment','left',...
                'verticalalignment','top');
        elseif j == max(ids)
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            h = text(xl(2),yl(2),sprintf('%d/%d',cells(k).p,cells(k).cl),...
                'horizontalalignment','right',...
                'verticalalignment','bottom');
        end
    end
    end
    end
    setappdata(DATA.toplevel,'xcorrs',xcorrs);
%    ReplotXcorrs(DATA, [], 'Shape/Efficacy')
    
function ReplotXcorrs(a,b, type)
       
    DATA = GetDataFromFig(a);
    xcorrs = getappdata(DATA.toplevel,'xcorrs');
    if sum(strcmpi(type, {'meanim' 'xcorrs' 'syncspikes' 'histograms'}))
        DATA.plot.xcorrtype = type;
        SetMenuCheck(a,'exclusive');
    elseif strcmpi(type, 'Shape/Efficacy')
        GetFigure(DATA.tag.xcorrpop);
        ClearPlot;
        hold off;
    for j = 1:length(xcorrs)
    plot(xcorrs(j).shapexc, max(xcorrs(j).efficacy),'o',...
        'buttondownfcn',{@PlotXcorr, xcorrs(j).p(1),xcorrs(j).p(2)});
    hold on;
    end
    
    set(gca,'yscale','log','ylim',[min([xcorrs.shapexc]) 1]);
    end

    
    function xc = xShapeCorr(P,Q)
    xc = corrcoef(P.MeanSpike.ms(:),Q.MeanSpike.ms(:));
    xc = xc(1,2);
