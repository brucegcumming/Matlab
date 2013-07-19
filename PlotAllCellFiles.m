function PlotAllCellFiles(name, varargin)

DATA.tag.allexpts = 'AllCellPlot';
DATA.toplevel = 0;

j = 1;
while j <= length(varargin)
    j = j+1;
end
if isstruct(name)
    DATA.Expts = name;
elseif isdir(name)
    exptlist = FindAllExpts(name);
    load(exptlist{1});
    DATA.exptlist = exptlist;
    DATA.Expts = AllExpt;
else
    load(name);
    DATA.Expts = AllExpt;
end
a = SetFigure(DATA.tag.allexpts,DATA);
DATA = get(a.toplevel,'UserData');
%set(DATA.toplevel,'UserData',DATA);
DATA = SetAllExpts(DATA);

function DATA = SetAllExpts(DATA)
DATA.AllExpts = {}; %clear any previous
for j = 1:length(DATA.Expts.Spikes)
    Expt = All2Expt(DATA.Expts,j,'all');
    DATA.AllExpts{j}.plotres = PlotExpt(Expt,'rcnmin',5);
    DATA.AllExpts{j}.Header = Expt.Header;
    DATA.AllExpts{j}.Stimvals = Expt.Stimvals;
    DATA.AllExpts{j}.cellid = Expt.Header.cellnumber;
end
PlotAllCells(DATA,[], 'rates');
setappdata(DATA.toplevel,'AllCellRes',DATA.AllExpts);

function PlotAllCells(a,b, plottype)
    
    DATA = GetDataFromFig(a);
    
   SetFigure(DATA.tag.allexpts,DATA,'force');
   peaks = [];
   hold off;
   if strncmp(plottype,'xcorr',5)
       PlotXcorrs(DATA, plottype);
       return;
   elseif strncmp(plottype,'offsetrateseq',10)
       PlotRateSequences(DATA,'offset','normalize');
       return;
   elseif strncmp(plottype,'normrateseq',10)
       PlotRateSequences(DATA,'normalize');
       return;
   elseif strncmp(plottype,'rateseq',7)
       PlotRateSequences(DATA);
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
       if isfield(Expts{j}.Header,'cellnumber')
           cells(j) = Expts{j}.Header.cellnumber;
       else
           cells(j) = 0;
       end
       if isfield(Expts{j}.plotres,'bestdelay')
           ny(j) = size(Expts{j}.plotres(1).y,2);
       elseif isfield(Expts{j}.plotres,'x')
           ny(j) = size(Expts{j}.plotres(1).x,2);
       end
       if isfield(Expts{j}.plotres,'x')
       allx = cat(1, allx, Expts{j}.plotres(1).x(:,1));
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
           P(1).Header.nspk = Expts{j}.Header.nspk;
           P(1).Header.dips = Expts{j}.Header.dips;
           if isfield(P,'Data')
               Expt = P(1).Data;
           end

%           id = find(ismember([Expts{j}.Trials.Trial],Expts{j}.Header.BlockStart));
%           P(1).Header.blocks = [Expts{j}.Trials(id).id];
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
           else
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
                       [a,tid] = ismember(P.ids{t},[Expt.Trials.id]);
                       if isfield(Expt.Trials,'xo')
                           P.xpos(t) = mean([Expt.Trials(tid).xo]);
                       end
                       if isfield(Expt.Trials,'yo')
                           P.ypos(t) = mean([Expt.Trials(tid).yo]);
                       end
                   end
               end
           end
           hold on;
           if isfield(P(1), 'fit') && isfield(P(1).fit,'xv')
               plot(P.fit.xv, P.fit.fitcurve.*scale,'color',colors{j});
               peaks(did(j)) = P.fit.peak;
           end
           AllCellRes(j,:) = P;
       end
   end
   mylegend(h,labels);

   X = DATA;
   X.sequence = did;
   set(gcf,'UserData',X);
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
           set(h,'buttondownfcn',{@HitImage, 'allexpt'});
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
   SetCellChooser(DATA);

   
function figpos = SetFigure(tag, DATA, varargin)
     
if DATA.toplevel > 0
    figpos = getappdata(DATA.toplevel,'FigPos');
else
    figpos.tag{1} = tag;
end

     j = 1;
     args = {};
     while j < length(varargin)
         args = {args{:} varargin{j}};
         j =  j+1;
     end
     [a,isnew] = GetFigure(tag,args{:});
     if isfield(figpos,'tag') && iscellstr({figpos.tag})
     id = strmatch(tag,{figpos.tag});
     if ~isempty(id) && (isnew || figpos(id).set ==0)
             set(a,'Position',figpos(id).pos);
             figpos(id).set = 1;
             setappdata(DATA.toplevel,'FigPos',figpos);
     end
     end
     if isnew && strcmp(DATA.tag.allexpts,tag)
         hm = uimenu(a,'label','&Plots');
         uimenu(hm,'label','&Next','Callback',{@AllCellPlots, 'next'});
         uimenu(hm,'label','&Prev','Callback',{@AllCellPlots, 'prev'});
         uimenu(hm,'label','Choose','Tag','ExptCellChoose');
         sm = uimenu(hm,'label','&Sort');
         uimenu(sm,'label','by &Var','Callback',{@AllCellPlots, 'sortbyvar'});
         uimenu(sm,'label','by &Probe','Callback',{@AllCellPlots, 'sortbyprobe'});
         uimenu(hm,'label','&Save AllExpts','Callback',{@AllCellPlots, 'save'});
         uimenu(hm,'label','&Xcorrs','Callback',{@AllCellPlots, 'xcorr'});
         set(a,'KeyPressFcn',@AllCellKeyPressed);
         if DATA.toplevel == 0
             DATA.toplevel = a;
             set(a,'UserData',DATA);
         else
             set(a,'UserData',DATA.toplevel);
         end
         figpos.toplevel = DATA.toplevel;
         if isfield(DATA,'exptlist')
             sm = uimenu(hm,'label','&Expts');
             for j = 1:length(DATA.exptlist)
                 uimenu(sm,'label',DATA.exptlist{j},'Callback',{@AllCellPlots, 'setexpt', j});
             end
         end
         SetCellChooser(DATA);
     elseif DATA.toplevel == 0
         DATA.toplevel = a;
         figpos.toplevel = a;
         set(a,'UserData',DATA);
     end

function SetCellChooser(DATA)
    it = findobj('Tag','ExptCellChoose');
    delete(get(it,'Children'));
    AllCellRes = getappdata(DATA.toplevel,'AllCellRes');
    
    for j = 1:length(AllCellRes)
        if AllCellRes(j).cellid > 0
            uimenu(it,'label',sprintf('Cell%d (%.1f)',AllCellRes(j).cellid,AllCellRes(j).Header.probe),'callback',{@AllCellPlots, j});
        else
            uimenu(it,'label',sprintf('P%dmu',AllCellRes(j).Header.probe),'callback',{@AllCellPlots, j});
        end
    end
    
function AllCellKeyPressed(src,ks, op)
     DATA = GetDataFromFig(src);
     AllCellRes = getappdata(DATA.toplevel,'AllCellRes');
     P =  getappdata(DATA.toplevel,'AllCellPlot');
     if strcmp(ks.Key,'downarrow') && P.currentcell < P.nplots
         P = NextPlot(P,AllCellRes,0);
     elseif strcmp(ks.Key,'uparrow') && P.currentcell > 1
         P = NextPlot(P,AllCellRes,-1);
     end
     setappdata(DATA.toplevel,'AllCellPlot',P);
     
function P =  NextPlot(P, AllCellRes, step)
    R = [];
    if P.sortbyvar
        var = cat(1,AllCellRes.var);
        id = find(var(:,2) < 0.4);
        var(id,2) = 0.4;
        [a, id] = sort(var(:,1)./var(:,2),'descend');
    else
        id = 1:P.nplots;
    end
    if step == 0 && P.currentcell < P.nplots
        hold off;
        P.currentcell = P.currentcell+1;
    elseif step == -1 && P.currentcell > 1
        hold off;
        P.currentcell = P.currentcell-1;
    else
        P.currentcell = step;
    end
    R = PlotResult(AllCellRes(id(P.currentcell)));
if P.sortbyvar && ~isempty(R)
    t = get(gca,'title');
    set(t,'string',[get(t,'string') sprintf('  VR %.3f (%.3f/%.3f)',R.var(1)./R.var(2),R.var(1),R.var(2))]);
end
     
function AllCellPlots(a,b, op, varargin)
     DATA = GetDataFromFig(a);
     j = 1;
     if strcmp(op,'setexpt')
         j = 2;
         load(DATA.exptlist{varargin{1}});
         DATA.Expts = AllExpt;
         DATA = SetAllExpts(DATA);
     end
     while j <= length(varargin)
         j = j+1;
     end
     onoff = {'off' 'on'};
     AllCellRes = getappdata(DATA.toplevel,'AllCellRes');
     if isfield(AllCellRes{1},'Header')
     ename = Expt2Name(AllCellRes{1});
     elseif isfield(DATA.Expts,'Spikes')
         ename = Expt2Name(DATA.Expts.Expt);
     else
         ename = Expt2Name(DATA.Expts{DATA.currentexpt});
     end
     P =  getappdata(DATA.toplevel,'AllCellPlot');
     if isnumeric(op)
         P = NextPlot(P,AllCellRes,op);
     elseif strcmp(op,'xcorr')
         F = gcf;
         oldname = get(F,'name');
         AllExpts = getappdata(DATA.toplevel,'AllExpt');
         set(F,'Name','Calculating Xcorrs');
         drawnow;
         xc = ExptListCorrs(AllExpts);
         ExptPlotXcorr(xc,1);
         set(F,'Name',oldname);
         drawnow;
     elseif strcmp(op,'save')
     outname = [DATA.Expt.Header.fileprefix '.Cellres.' ename '.mat'];
     save(outname,'AllCellRes');
     fprintf('Saved Cell Results to %s\n',outname');
     elseif sum(strcmp(op,{'sortbyvar' 'sortbyprobe'}))
         P.(op) = ~P.(op);
         set(a,'checked',onoff{1+P.(op)});
         
     elseif strcmp(op,'next')
         P = NextPlot(P,AllCellRes,0);
     elseif strcmp(op,'prev')
         P = NextPlot(P,AllCellRes,-1);
     end
     setappdata(DATA.toplevel,'AllCellPlot',P);

    
function exptlist = FindAllExpts(name)

d = dir([name '/*Cells*.mat']);
for j = 1:length(d)
    exptlist{j} = [name '/' d(j).name];
end

        