function PlotAllCellFiles(name, varargin)
% plot expts from an allexpt Structure

% makes appdata AllCellRes, with just plot results
% Expts with cell array of  of Expt structs where possible
DATA.tag.allexpts = 'AllCellPlot';
DATA.toplevel = 0;
plotargs = {};
plottype = 'rates';
parentfigure = 0;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'efficacy',5)
        plottype = varargin{j};
    elseif strncmpi(varargin{j},'parent',5)
        j = j+1;
        parentfigure = varargin{j};
    elseif strncmpi(varargin{j},'tag',5)
        j = j+1;
        DATA.tag.allexpts = varargin{j};
    end
    plotargs{j} = varargin{j};
    j = j+1;
end

if isstruct(name) %AllExpt File
    AllExpt = expt.CheckFields(name);
elseif iscell(name) %online AllExpt Struct
    DATA.AllExpts = CheckAllRes(name);
elseif isdir(name)
    exptlist = FindAllExpts(name);
    load(exptlist{1});
    DATA.exptlist = exptlist;
else
    load(name); %?load Expt?
end
[a, F] = SetFigure(DATA.tag.allexpts,DATA);
if double(parentfigure) > 0
    setappdata(F, 'ParentFigure', parentfigure);
end

DATA = get(a.toplevel,'UserData');
if ~isfield(DATA,'AllExpts')
    
    DATA.Expts = CheckAllExpt(AllExpt);
%    setappdata(DATA.toplevel,'Expts',AllExpt);
    
    %set(DATA.toplevel,'UserData',DATA);
    if strncmp(plottype,'effic',5)
        PlotAllEfficacies(DATA);
    else
        DATA = SetAllExpts(DATA, plotargs{:});
    end
else
    PlotAllCells(DATA,[], 'rates');
    setappdata(F,'Expts',DATA.AllExpts);
    if isfigure(parentfigure)%not sure we need this.
        %setappdata(parentfigure,'AllCellRes',DATA.AllExpts);
    end
end

function A = CheckAllRes(A)

for j = 1:length(A)
    if isfield(A{j},'Header') && isfield(A{j}.plotres,'Header')
        A{j}.plotres(1).Header = CopyFields(A{j}.plotres(1).Header,A{j}.Header,'probe','cellnumber');
    end
end


function A = CheckAllExpt(A)

for j = 1:length(A.Header)
    if ~isfield(A.Header,'quality')
        A.Header(j).quality = 0;
    end
end

function DATA = SetAllExpts(DATA, varargin)
AllExpts = {}; %clear any previous
args = {};
parfor j = 1:length(DATA.Expts.Spikes)
    Expt = All2Expt(DATA.Expts,j,'all');
    AllExpts{j}.plotres = PlotExpt(Expt,'rcnmin',5,'fbox','noplot',varargin{:});
    AllExpts{j}.Header = Expt.Header;
    AllExpts{j}.Trials = Expt.Trials;
    AllExpts{j}.Stimvals = Expt.Stimvals;
    AllExpts{j}.cellid = Expt.Header.cellnumber;
    AllExpts{j}.plotres(1).cellid = Expt.Header.cellnumber;
end
DATA.AllExpts = AllExpts;
PlotAllCells(DATA,[], 'rates');
setappdata(DATA.toplevel,'Expts',DATA.AllExpts);

function PlotAllEfficacies(DATA, varargin)
for j = 1:length(DATA.Expts.Spikes)
    Expts{j} = All2Expt(DATA.Expts,j,'all');
    spkt{j} = GetSpkTimes(Expts{j});
    labels{j} = sprintf('%.1f',DATA.Expts.Header(j).probe);
end
setappdata(gcf,'SpkTimes',spkt);
DATA.efficacies = CalcEfficacies(spkt);
SetData(DATA);
hold off;
imagesc(DATA.efficacies,'buttondownfcn',@HitImage);
set(gca,'ytick',1:length(labels),'yticklabel',labels);
set(gca,'xtick',1:length(labels),'xticklabel',labels);
colorbar;


function PlotAllXCorr(DATA, varargin)
for j = 1:length(DATA.Expts.Spikes)
    Expts{j} = All2Expt(DATA.Expts,j,'all');
    spkt{j} = GetSpkTimes(Expts{j});
    labels{j} = sprintf('%.1f',DATA.Expts.Header(j).probe);
end
setappdata(gcf,'SpkTimes',spkt);
DATA.efficacies = CalcEfficacies(spkt);
SetData(DATA);
hold off;
imagesc(DATA.efficacies,'buttondownfcn',@HitImage);
set(gca,'ytick',1:length(labels),'yticklabel',labels);
set(gca,'xtick',1:length(labels),'xticklabel',labels);
colorbar;


function HitLine(a,b, varargin)

p = varargin{1};
    DATA = get(gcf,'UserData');
    X = getappdata(gcf,'PlotData');
    A = getappdata(gcf,'AllCellRes');

    if isfield(DATA,'Expts')
        PlotExpt(All2Expt(DATA.Expts,p,'all'));
    else
        PlotResult(A(p));
    end

function HitImage(a,b, varargin)

pos = get(gca,'currentpoint');
x = round(pos(1,[1 2]));
plottype = 'ccf';
j = 1;
while j <= length(varargin)
    if strcmp(varargin{j},'allexpt')
        plottype = 'expt';
    elseif strcmp(varargin{j},'rateseq')
        plottype = 'rateseq';
    elseif strcmp(varargin{j},'xcorr')
        plottype = 'xcorr';
    end
    j = j+1;
end
    DATA = get(gcf,'UserData');
    X = getappdata(gcf,'PlotData');
    A = getappdata(gcf,'AllCellRes');
    AllExpt = getappdata(gcf,'AllExpts');
    Expts = getappdata(gcf,'Expts');
    if isempty(Expts)
        Expts = DATA.AllExpts;
    end
    eid = x(2);
    if length(Expts) >= X.depthorder(eid)
        E = Expts{X.depthorder(eid)};
    end
    
if strcmp(plottype,'expt')
    GetFigure('Data');
    eid = x(2);
    if isfield(DATA,'Expts') && isfield(DATA.Expts,'Expt') %AllExpt Struct
        PlotExpt(All2Expt(DATA.Expts,X.depthorder(eid),'all'));
    elseif ~isempty(Expts)
        PlotExpt(E);
    else
        PlotResult(A(X.depthorder(eid)),'ParentFigure',DATA.toplevel);
    end
elseif strcmp(plottype,'rateseq')
    GetFigure('Data');
    hold off; 
    PlotRateSequence(E);
    if E.Header.cellnumber > 0
    else
        title(sprintf('P%d',E.Header.probe));
    end
elseif strcmp(plottype,'xcorra')
else
    spkt = getappdata(gcf,'SpkTimes');
    if ~isempty(Expts)
        A = Expts{x(2)}.Header;
        B = Expts{x(1)}.Header;
        [xc, details] = expt.xcorr(Expts([x(2) x(1)]));
        t = details.xpts;
        str = sprintf('P%d->%d',A.probe,B.probe);
    elseif ~isempty(AllExpts)
        A = AllExpts.Header(x(2));
        B = AllExpts.Header(x(1));
    str = sprintf('C%d/P%dQ%.1f vs C%d/P%dQ%.1f,',A.cellnumber,A.probe,A.quality,B.cellnumber,B.probe,B.quality);
    fprintf('Efficacy %d-%d %.3f (%.3f)\n',A.probe,B.probe,DATA.efficacies(x(2),x(1)),DATA.efficacies(x(1),x(2)));
    fprintf('%s\n',str);
    GetFigure('xcorr');
    [xc, details] = xcorrtimes(spkt{x(2)}./10000,spkt{x(1)}./10000);
    t = details.t;
    end
    GetFigure('xcorr');
    plot(t,xc);
    title(str);
end


function PlotAllCells(a,b, plottype)
     
   if ishandle(a)
       DATA = get(GetFigure(a),'UserData');
   else
       DATA = a;
   end
    
   [~,F] = SetFigure(DATA.tag.allexpts,DATA,'force');
   peaks = [];
   hold off;
   if strncmp(plottype,'xcorr',5)
       PlotXcorrs(DATA, plottype);
       return;
   elseif strncmp(plottype,'offsetrateseq',10)
       PlotRateSequences(DATA,plottype);
       return;
   elseif strncmp(plottype,'normrateoffset',12)
       PlotRateSeqIm(DATA);
       PlotRateSequences(DATA,plottype);
       return;
   elseif strncmp(plottype,'normrateseq',12)
       PlotRateSeqImage(DATA);
%       PlotRateSequences(DATA,plottype);
       return;
   elseif strncmp(plottype,'rateseq',7)
       PlotRateSequences(DATA,plottype);
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
   ytype = '';
   ztype = '';
   allx = [];
   ally = [];
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
       if isfield(Expts{j}.plotres,'y')
           ally = cat(1, ally, Expts{j}.plotres(1).y(1,:)');
       end
   end
   collapse = [0 0 1];
   if sum(strcmp(Expts{1}.Stimvals.e3,{'ce' 'a2' 'mixac'}))
       collapse = [0 0 0];
   end
   allx = unique(allx);
   ally = unique(ally);
   Aplot = getappdata(DATA.toplevel,'AllCellPlot');
   if isfield(Aplot,'sortcellsbyprobe') && Aplot.sortcellsbyprobe %only cells
       id = find(cells > 0);
       depths = depths(id);
       Expts = Expts(id);
       ny = ny(id);
       cells = cells(id);
   end

   ny = max(ny);
   [a,did] = sort(depths);
   X.depthorder = did;
   setappdata(gcf,'PlotData',X);
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
           if isfield(Expts{j}.Header,'dips')               
           P(1).Header.dips = Expts{j}.Header.dips;
           end
           if isfield(P,'Data')
               Expt = P(1).Data;
           end

%           id = find(ismember([Expts{j}.Trials.Trial],Expts{j}.Header.BlockStart));
%           P(1).Header.blocks = [Expts{j}.Trials(id).id];
           if strmatch(plottype,{'sacsdf' 'sacsdfim'})
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
                   h = plot(V.*scale ,'-','color',colors{j});
                   im(row,1:length(V),yi) = V.*scale;
               elseif strmatch(plottype,{'blank' 'blankim'})
                   xi = find(P.sdfs.extraval == -1009);
                   V = P.sdfs.extras{xi}.sdf(ts:end);
                   scale = 1./median(V);
                   h = plot(V.*scale,'-','color',colors{j});
                   im(row,1:length(V),yi) = V.*scale;
               else
                   h = plot(P.x(:,t),P.y(:,yi,t).*scale,'o-','color',colors{j});
                   im(row,iid,1) = Expts{j}.plotres.y(:,1,t).*scale;
               end
               set(h,'buttondownfcn',{@HitLine, j});
               hold on;
           elseif isfield(P(1),'means')
               P(1).var(1) = std(sqrt(P(1).means(:)));
               vars = [];
               for k = 1:length(P(1).counts(:))
                   vars(k) = var(sqrt(P(1).counts{k}));
               end
               P(1).var(2) = sqrt(mean(vars));
               if collapse(3)
                   P(1).x = mean(P(1).x,3);
                   P(1).y = mean(P(1).y,3);
                   P(1).n = sum(P(1).n,3);
                   P(1).means = mean(P(1).means,3);
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
               hold on;
           else
               fprintf('Missing Data/No Spikes for P%d\n',j);
           end
           if isfield(P(1),'type') && length(P(1).type) > 1 && ~isempty(P(1).type{2})
               ytype = P(1).type{2};
           end
           if isfield(P(1),'type') && length(P(1).type) > 2 && ~isempty(P(1).type{3})
               ztype = P(1).type{3};
           end
           if isfield(P(1), 'fit') && isfield(P(1).fit,'xv')
               plot(P.fit.xv, P.fit.fitcurve.*scale,'color',colors{j});
               peaks(did(j)) = P.fit.peak;
           end
           try
           AllCellRes(j,:) = P;
           catch
               for k = 1:length(P)
                   f = fields(P(k));
                   for c = 1:length(f)
                       AllCellRes(j,k).(f{c}) = P(k).(f{c});
                   end                   
               end
           end
       end
   end
   mylegend(h,labels);

   X = DATA;
   X.sequence = did;
   set(gcf,'UserData',X);
   %set defaults then copy in existing values
   AllCellPlot.nplots = length(AllCellRes);
   AllCellPlot.currentcell = 0;
   AllCellPlot.sortbyvar = 0;
   AllCellPlot.sortbyprobe = 1;
   AllCellPlot.sortcellsbyprobe = 0;
   AllCellPlot = CopyFields(AllCellPlot,Aplot);
   setappdata(F,'AllCellPlot',AllCellPlot);
   setappdata(F,'AllCellRes',AllCellRes);
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
           title(sprintf('%s=%.2f',ytype,ally(j)));
           if sum(cells ==0) == 0
               set(gca,'ytick',[]);
               for k = 1:size(im,1)
                   text(xl(1),k,sprintf('%d:%.1f',Expts{did(k)}.Header.cellnumber,depths(did(k))));
               end
           end
       end
   end
   SetCellChooser(F);

   
function [figpos, F] = SetFigure(tag, DATA, varargin)
     
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
     [F,isnew] = GetFigure(tag,args{:});
     if isfield(figpos,'tag') && iscellstr({figpos.tag})
     id = strmatch(tag,{figpos.tag});
     if ~isempty(id) && (isnew || figpos(id).set ==0)
             set(F,'Position',figpos(id).pos);
             figpos(id).set = 1;
             setappdata(DATA.toplevel,'FigPos',figpos);
     end
     end
     if isnew && strcmp(DATA.tag.allexpts,tag)

         hm = uimenu(F,'label','&Plots','Tag','PlotMenu');
         uimenu(hm,'label','&Next','Callback',{@AllCellPlots, 'next'});
         uimenu(hm,'label','&Prev','Callback',{@AllCellPlots, 'prev'});
         uimenu(hm,'label','Choose','Tag','ExptCellChoose');
         sm = uimenu(hm,'label','&Type');
         uimenu(hm,'label','Rate &Seq','Callback',{@PlotAllCells, 'rateseq'});
         uimenu(hm,'label','Rates','Callback',{@PlotAllCells, 'rates'});
         uimenu(hm,'label','Normalized (Max) Rates','Callback',{@PlotAllCells, 'normmax'});
         uimenu(hm,'label','Normalized (Max) Rate image','Callback',{@PlotAllCells, 'imagemax'});
         uimenu(hm,'label','Normalized (Mean) Rates','Callback',{@PlotAllCells, 'normmean'});
         uimenu(hm,'label','Normalized (Mean) Rate image','Callback',{@PlotAllCells, 'imagemean'});
         uimenu(hm,'label','Norm Rate &Seq','Callback',{@PlotAllCells, 'normrateseq'});
         uimenu(hm,'label','Rate &Seq+offset','Callback',{@PlotAllCells, 'offsetrateseq'});
         uimenu(hm,'label','Test','Callback',{@PlotAllCells, 'test'});
         sm = uimenu(hm,'label','&Sort');
         uimenu(sm,'label','by &Var','Callback',{@AllCellPlots, 'sortbyvar'});
         uimenu(sm,'label','by &Probe','Callback',{@AllCellPlots, 'sortbyprobe'});
         uimenu(sm,'label','Cells only by &Probe','Callback',{@AllCellPlots, 'sortcellsbyprobe'});
         uimenu(hm,'label','&Save AllExpts','Callback',{@AllCellPlots, 'save'});
         uimenu(hm,'label','&Xcorr efficacy','Callback',{@AllCellPlots, 'efficacy'});
         uimenu(hm,'label','&Xcorrs','Callback',{@AllCellPlots, 'xcorr'});
         set(F,'KeyPressFcn',@AllCellKeyPressed);
         if DATA.toplevel == 0
             DATA.toplevel = F;
             set(F,'UserData',DATA);
         else
             set(F,'UserData',DATA.toplevel);
         end
         figpos.toplevel = DATA.toplevel;
         if isfield(DATA,'exptlist')
             sm = uimenu(hm,'label','&Expts');
             for j = 1:length(DATA.exptlist)
                 uimenu(sm,'label',DATA.exptlist{j},'Callback',{@AllCellPlots, 'setexpt', j});
             end
         end
         if isfield(DATA,'AllExpts')
%             setappdata(F,'AllCellRes',DATA.AllExpts);
%             SetCellChooser(F);
         end
     elseif DATA.toplevel == 0
         DATA.toplevel = F;
         figpos.toplevel = F;
         set(F,'UserData',DATA);
     end

function SetCellChooser(F)
    it = findobj(F,'Tag','ExptCellChoose'); 
    if ~isempty(it)
        delete(get(it,'Children'));
    else
        it = findobj(allchild(F),'flat','Tag','PlotMenu');
        it = uimenu(it,'label','Choose','Tag','ExptCellChoose');
    end
    AllCellRes = getappdata(F,'AllCellRes');
    
    for j = 1:length(AllCellRes)
        if AllCellRes(j).cellid > 0
            uimenu(it,'label',sprintf('Cell%d (%.1f)',AllCellRes(j).cellid,AllCellRes(j).Header.probe),'callback',{@AllCellPlots, j});
        else
            uimenu(it,'label',sprintf('P%dmu',AllCellRes(j).Header.probe),'callback',{@AllCellPlots, j});
        end
    end
    
function AllCellKeyPressed(src,ks, op)

     F = GetFigure(src);
     DATA = GetDataFromFig(src);
     AllCellRes = getappdata(F,'AllCellRes');
     P =  getappdata(F,'AllCellPlot');
     if strcmp(ks.Key,'downarrow') && P.currentcell < P.nplots
         P = NextPlot(P,AllCellRes,0);
     elseif strcmp(ks.Key,'uparrow') && P.currentcell > 1
         P = NextPlot(P,AllCellRes,-1);
     end
     setappdata(F,'AllCellPlot',P);
set(0,'CurrentFigure',gcbf);
     
     
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
    hold off;
    R = PlotResult(AllCellRes(id(P.currentcell)));
    t = get(gca,'title');
    xstr = '';
    if isfield(R(1).Header,'probe')
        xstr = sprintf('P%d',R(1).Header.probe);
    end

if P.sortbyvar && ~isempty(R)
    set(t,'string',[get(t,'string') sprintf('  VR %.3f (%.3f/%.3f)',R.var(1)./R.var(2),R.var(1),R.var(2))]);
else
    set(t,'string',[get(t,'string') xstr]);
end
     
function AllCellPlots(a,b, op, varargin)
     F = GetFigure(a);
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
     AllCellRes = getappdata(F,'AllCellRes');
     Expts = getappdata(F,'Expts');
     AllExpts = getappdata(F,'AllExpt');
     Header = [];
     if iscell(AllCellRes)
         for j = 1:length(AllCellRes)
             if isempty(Header) && isfield(AllCellRes{j},'Header') 
                 Header = AllCellRes{j}.Header;
             end
         end
     elseif isstruct(AllCellRes) && isfield(AllCellRes,'Header')
         for j = 1:length(AllCellRes)
             if isempty(Header) && isfield(AllCellRes(j).Header,'Name')
                 Header = AllCellRes(j).Header;
             end
         end
     end
     if ~isempty(Header)
         ename = Expt2Name(Header);
     elseif isfield(DATA.Expts,'Spikes')
         ename = Expt2Name(DATA.Expts.Expt);
     else
         ename = Expt2Name(DATA.Expts{DATA.currentexpt});
     end
     P =  getappdata(F,'AllCellPlot');
     if isnumeric(op)
         P = NextPlot(P,AllCellRes,op);
     elseif strcmp(op,'efficacy')
         if ~isempty(Expts)
             I = CalcEfficacies(Expts);
             hold off;
             h = imagesc(I);
             xlabel('Probe A');
             ylabel('Probe B');
             title('P(simultaneous');
             colorbar;
             set(h,'buttondownfcn',{@HitImage, 'xcorr'});
         end
     elseif strcmp(op,'xcorr')
         F = gcf;
         oldname = get(F,'name');
         AllExpts = getappdata(F,'AllExpt');
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
     elseif sum(strcmp(op,{'sortbyvar' 'sortbyprobe' 'sortcellsbyprobe'}))
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


function PlotRateSeqImage(DATA, type, varargin)
Expts = getappdata(DATA.toplevel,'Expts');
offset = 0;
hold off;
alltimes = [];
for j = 1:length(Expts)
    X = PlotRateSequence(Expts{j},'normalize','noplot');
    rateim(j,X.times) = X.rates;
    alltimes = [alltimes X.times];
end
altimes = unique(alltimes);
h = imagesc(rateim);
set(h,'buttondownfcn',{@HitImage, 'rateseq'});
xlabel('trial');
ylabel('probe');
title('trial rates');

function PlotRateSequences(DATA, type, varargin)

colors = mycolors;
Expts = getappdata(DATA.toplevel,'Expts');
offset = 0;
hold off;
for j = 1:length(Expts)
    if strcmp(type,'offsetrateseq')
        args = {'offset' offset};
    elseif strcmp(type,'normrateseq')
        args = {'normalize' 'offset' j};
    else
        args = {};
    end
    X = PlotRateSequence(Expts{j},args{:},'color',colors{j});
    offset = offset+max(X.rates);
    hold on;
end