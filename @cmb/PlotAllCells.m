function PlotAllCells(a,b, plottype)

DATA = GetDataFromFig(a);

cmb.SetFigure(DATA.tag.allexpts,DATA,'force');
peaks = [];
hold off;
if strncmp(plottype,'xcorr',5)
    cmb.PlotXcorrs(DATA, plottype);
    return;
elseif strncmp(plottype,'offsetrateseq',10)
    cmb.PlotRateSequences(DATA,'offset','normalize');
    return;
elseif strncmp(plottype,'normrateseq',10)
    cmb.PlotRateSequences(DATA,'normalize');
    return;
elseif strncmp(plottype,'rateseq',7)
    cmb.PlotRateSequences(DATA);
    return;
elseif strncmp(plottype,'load',4)
    outname = [DATA.Expt.Header.fileprefix '.Cellres.' DATA.Expt.Header.expname '.mat'];
    load(outname);
    AllCellPlot.nplots = length(AllCellRes);
    AllCellPlot.currentcell = 0;
    setappdata(DATA.toplevel,'AllCellPlot',AllCellPlot);
    setappdata(DATA.toplevel,'AllCellRes',AllCellRes);
    return;
end
Expts = DATA.AllExpts;
colors = mycolors;
lines = {'-','--',':','-.','-','--',':','-.','-','--',':','-.','-','--',':','-.'};
h = []; labels = {};
allx = [];
for j = 1:length(Expts)
    depths(j) = Expts{j}.Header.probe;
    if ~isfield(Expts{j}.Header,'cellnumber')
        Expts{j}.Header.cellnumber = 0;
    end
    cells(j) = Expts{j}.Header.cellnumber;
    good = 1;
    if isfield(Expts{j}.plotres,'bestdelay')
        ny(j) = size(Expts{j}.plotres(1).y,2);
    elseif isfield(Expts{j}.plotres,'x')
        ny(j) = size(Expts{j}.plotres(1).x,2);
    else
        good = 0;
    end
    if good
        allx = cat(1, allx, Expts{j}.plotres(1).x(:,1));
    else
    end
end
collapse = [0 0 1];
if sum(strcmp(Expts{1}.Stimvals.e3,{'ce' 'a2' 'mixac'}))
    collapse = [0 0 0];
end
allx = unique(allx);
ny = max(ny);
[a,did] = sort(depths);
ts = 20; %ignore first 20 samples of sdf
normtype = 1;
subplot(1,1,1);
allxv = [];
for yi = 1:ny;
    if strmatch(plottype,{'normmax' 'imagemax'})
        im(1:length(Expts),1:length(allx),yi) = NaN;
    end
    for j = 1:length(Expts)
        row = find(did == j);
        P = Expts{j}.plotres;
        P(1).cellid = cells(j);
        P(1).Header.probe = Expts{j}.Header.probe;
        P(1).Header = CopyFields(P(1).Header,Expts{j}.Header,'nspk', 'dips');
        
        id = find(ismember([Expts{j}.Trials.Trial],Expts{j}.Header.BlockStart));
        P(1).Header.blocks = [Expts{j}.Trials(id).id];
        if strmatch(plottype,{'test' 'sacsdf' 'sacsdfim'})
            sdf = mksacsdf(Expts{j}.Trials,100,[-500:10:2000],'minsize',1);
            plot(sdf,'color',colors{j});
            if 0
                im(row,:) = sdf;
            else
                im(row,:) = sdf./mean(sdf(1:30));
            end
        elseif isfield(P,'bestdelay')
            scale = 1;
            if strmatch(plottype,{'normmax' 'imagemax'})
                scale = 1./max(P.y(:,P.bestdelay));
            elseif strmatch(plottype,{'normmean' 'imagemean'})
                scale = 1./mean(P.y(:,P.bestdelay));
            else
                scale = 1;
            end
            t = P.bestdelay;
            allxv = cat(1,allxv,  P.x(:,1));
            iid = find(ismember(allx, P.x(:,1)));
            if strmatch(plottype,{'var' 'varim'})
                V = std(cat(2,P.sdfs.s{:}),[],2);
                if normtype == 1
                    scale = 1./median(V);
                end
                plot(V.*scale ,'-','color',colors{j});
                im(row,1:length(V),yi) = V.*scale;
            elseif strmatch(plottype,{'blank' 'blankim'})
                xi = find(P.sdfs.extraval == -1009);
                V = P.sdfs.extras{xi}.sdf(ts:end);
                scale = 1./median(V);
                plot(V.*scale,'-','color',colors{j});
                im(row,1:length(V),yi) = V.*scale;
            else
                plot(P.x(:,t),P.y(:,yi,t).*scale,'o-','color',colors{j});
                im(row,iid,1) = Expts{j}.plotres.y(:,1,t).*scale;
            end
        elseif isfield(P(1),'means')
            P(1).var(1) = std(sqrt(P(1).means(:)));
            vars = [];
            for k = 1:length(P(1).counts(:))
                vars(k) = var(sqrt(P(1).counts{k}));
            end
            P(1).var(2) = sqrt(mean(vars));
            if collapse(3)
                P(1).x = sum(P(1).x,3);
                P(1).y = sum(P(1).y,3);
                P(1).n = sum(P(1).n,3);
                P(1).means = sum(P(1).means,3);
            end
            nmin = 1+ median(P(1).n(:))/3;
            aid = find(P(1).n(:) >= nmin);
            if size(P(1).n,2) >= yi %some cells might be missing some values
                [xid, zid] = find(P(1).n(:,yi,:) >= nmin);
                if strmatch(plottype,{'normmax' 'imagemax'})
                    scale = 1./max(P(1).means(aid));
                elseif strmatch(plottype,{'normmean' 'imagemean'})
                    scale = 1./mean(P(1).means(aid));
                else
                    scale = 1;
                end
                xv = P(1).x(xid,yi,1);
                for zi = unique(zid)'
                    h(j) = plot(xv,squeeze(P(1).means(xid,yi,zi)).*scale,'o-','color',colors{j},'linestyle',lines{yi});
                end
                labels{j} = sprintf('C%d',Expts{j}.Header.cellnumber);
                if strmatch(plottype,{'imagemax' 'imagemean'  })
                    iid = find(ismember(allx, xv));
                    im(row,iid,yi) = P(1).means(xid,yi,1).*scale;
                end
            end
            allxv = [allxv;  P(1).x(:,1)];
            if strmatch(P(1).type{1},{'Pp' 'Op'})
                for t = 1:length(P.ids);
                    [a,tid] = ismember(P.ids{t},[Expts{j}.Trials.id]);
                    if isfield(Expts{j}.Trials,'xo')
                        P.xpos(t) = mean([Expts{j}.Trials(tid).xo]);
                    end
                    if isfield(Expts{j}.Trials,'yo')
                        P.ypos(t) = mean([Expts{j}.Trials(tid).yo]);
                    end
                end
            end
        else
            fprintf('Missing data for Expt %d\n',j);
        end
        hold on;
        if isfield(P(1), 'fit') && isfield(P(1).fit,'xv')
            plot(P.fit.xv, P.fit.fitcurve.*scale,'color',colors{j});
            peaks(did(j)) = P.fit.peak;
        end
        try
            AllCellRes(j,:) = P;
        catch
            pf = fields(P);
            af = fields(AllCellRes);
            xf = setdiff(pf,af);
            if ~isempty(xf)
                fprintf('%s new fields%s\n',IDString(P),sprintcell(xf,' %s'));
            end
            xf = setdiff(af,pf);
            if ~isempty(xf)
                fprintf('%s Missing fields%s\n',IDString(P),sprintcell(xf,' %s'));
            end
            for k = 1:size(P)
                AllCellRes(j,k).name = P(k).name;
                f = fields(P(k));
                for c = 1:length(f)
                    AllCellRes(j,k).(f{c}) = P(k).(f{c});
                end
            end
        end
    end
end
mylegend(h,labels);

X.toplevel = DATA.toplevel;
X.sequence = did;
set(gcf,'UserDATA',X);
AllCellPlot.nplots = length(AllCellRes);
AllCellPlot.currentcell = 0;
AllCellPlot.sortbyvar = 0;
AllCellPlot.sortbyprobe = 1;
setappdata(DATA.toplevel,'AllCellPlot',AllCellPlot);
setappdata(DATA.toplevel,'AllCellRes',AllCellRes);
xv = unique(allxv(:));
if strmatch(plottype,{'imagemax' 'imagemean' 'varim' 'blankim' 'test' 'sacsdfim'},'exact')
    for j =  1:ny
        subplot(1,ny,j);
        
        hold off;
        h = imagesc(minmax(xv),1:length(depths),squeeze(im(:,:,j)));
        set(h,'buttondownfcn',{@cmb.HitImage, 'allexpt'});
        caxis(minmax(im(:)));
        xl = get(gca,'xlim');
        if length(peaks)
            hold on;
            plot(peaks,1:length(peaks),'w');
        end
        if sum(cells ==0) == 0
            set(gca,'ytick',[]);
            for k = 1:size(im,1)
                text(xl(1),k,sprintf('%d:%.1f',Expts{did(k)}.Header.cellnumber,depths(did(k))));
            end
        end
    end
end
cmb.SetCellChooser(DATA);

