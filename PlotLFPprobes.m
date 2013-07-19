function result = PlotLFPprobes(rc, varargin)
%plot an LPF result stucture from PlotREvCOrAny 
%for multiple electrode recodings.
plottype = 1;
mtimes = [];
result = [];
xid = 1;
addblank = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'lines',5)
        plottype = 2;
    elseif strncmpi(varargin{j},'+blank',5)
        addblank =1;
    elseif strncmpi(varargin{j},'images',5)
        plottype = 4;
    elseif strncmpi(varargin{j},'eigscatter',5)
        plottype = 6;
    elseif strncmpi(varargin{j},'eigresp',5)
        plottype = 3;
    elseif strncmpi(varargin{j},'eigbresp',5)
        plottype = 9;
    elseif strncmpi(varargin{j},'eigb',4)
        plottype = 7;
    elseif strncmpi(varargin{j},'eiga',4)
        plottype = 8;
    elseif strncmpi(varargin{j},'onestim',5)
        plottype = 5;
    elseif strncmpi(varargin{j},'times',5)
        j = j+1;
        mtimes = varargin{j};
    elseif strncmpi(varargin{j},'xid',3)
        j = j+1;
        xid = varargin{j};
    end
    j = j+1;
end

if isfield(rc.sdfs,'z') % second order
    if length(xid) > 1
        for j = 1:length(xid)
            lfpresp(:,:,:,j) = cell2mat(rc.sdfs.lfp(xid(j),:,:));
        end
        res.lfpresp = mean(lfpresp,4);
    else
   res.lfpresp = cell2mat(rc.sdfs.lfp(xid,:,:));
    end
   res.lfpresp = res.lfpresp-repmat(rc.sdfs.alllfp,[1 1 size(res.lfpresp,3)]);
   if ~isempty(mtimes)
       id = find(rc.sdfs.lfptimes > mtimes(1) & rc.sdfs.lfptimes < mtimes(2));
       res.lfpresp = res.lfpresp(id,:,:);
       res.lfptimes = rc.sdfs.lfptimes(id);
   else
       id = find(rc.sdfs.lfptimes > 0);
       res.lfptimes = rc.sdfs.lfptimes(id);
   end
else
   res.lfpresp = cat(3,rc.sdfs.lfp{:});
   res.lfpresp = res.lfpresp-repmat(rc.sdfs.alllfp,[1 1 size(res.lfpresp,3)]);
   if ~isempty(rc.sdfs.extras)
       
       res.lfpblank = rc.sdfs.extras{1}.lfp-rc.sdfs.alllfp;
       if addblank
       res.lfpresp = cat(3, res.lfpblank,res.lfpresp);
       end
   end
   if ~isempty(mtimes)
       id = find(rc.sdfs.lfptimes > mtimes(1) & rc.sdfs.lfptimes < mtimes(2));
       res.lfpresp = res.lfpresp(id,:,:);
       if isfield(res,'lfpblank')
           res.lfpblank = res.lfpblank(id,:);
       end
       res.lfptimes = rc.sdfs.lfptimes(id);
   else
       id = find(rc.sdfs.lfptimes > 0);
       res.lfptimes = rc.sdfs.lfptimes(id);
   end
end
   res.lineplot = 0;
   res.plot = 0;
   res.x = squeeze(rc.sdfs.x(1,:,:));
   res.probes = 1:8;
   if plottype == 1
   PlotMovie(res, rc, 0, rc.sdfs.lfptimes(id));
   elseif plottype == 2
   PlotByProbe(res, rc,0);
   elseif ismember(plottype,[3 6 7 8 9])
   [result.resps,b] = PlotProbeEig(res, rc, plottype);
   result.blankresps = b.blankresps;
   elseif plottype == 4
   PlotByProbe(res, rc,1);
   elseif plottype == 5
   PlotByStim(res, rc,1);
   end


function [resps, details] = PlotProbeEig(res, rc, type, varargin)


lfpblank = [];
j = 1;
while j <= length(varargin)
    if strnmcpi(varargin{j},'blank')
        j = j+1;
        lfpblank = varargin{j};
    end
    j =j+1;
end

for probe = 1:length(res.probes)
    subplot(2,4,probe);
    [A,B] = eig(cov(squeeze(res.lfpresp(:,probe,:))'));
    for j = 1:10
        alleigs(j,:,probe) = A(:,end-j+1);
    end
    eigs(:,probe) = A(:,end);
    beigs(:,probe) = A(:,end-1);
    ceigs(:,probe) = A(:,end-2);

end
    if isfield(res,'lfpblank')
        for j = 1:size(alleigs,1)
            blankresps(j,:) = sum(res.lfpblank .* squeeze(alleigs(j,:,:)));
        end
        for j = 1:size(blankresps,2)
            for k = 1:size(blankresps,1);
            alleigs(k,:,j) = alleigs(k,:,j) .* sign(blankresps(k,j));
            end
            eigs(:,j) = eigs(:,j) .* sign(blankresps(1,j));
            beigs(:,j) = beigs(:,j) .* sign(blankresps(2,j));
        end
        blankresps = blankresps .* sign(blankresps);
    end
    for k = 1:size(res.lfpresp,3)
    xresps(k,:) = sum(res.lfpresp(:,:,k) .* eigs);
    yresps(k,:) = sum(res.lfpresp(:,:,k) .* beigs);
    for j = 1:size(alleigs,1)
        allresps(k,:,j) = sum(squeeze(res.lfpresp(:,:,k)) .* squeeze(alleigs(j,:,:)));
    end
    end
 subplot(1,1,1);
 hold off;
 if type == 6
 colors = mycolors;
 for j = 1:size(xresps,1)
     plot(xresps(j,:),yresps(j,:),'o','color',colors{j});
     hold on;
 end
 elseif type == 7
 imagesc(yresps);
 elseif type == 8
 imagesc(xresps);
 elseif type == 9
 imagesc(beigs);
 else
 imagesc(eigs);
 end
 
 resps = cat(3,xresps,yresps);
 resps = allresps;
 details.blankresps = blankresps;
 
function PlotByProbe(res, rc, type)

for j = 1:length(res.probes)
    subplot(2,4,j);
    if type == 1
    imagesc(squeeze(res.lfpresp(:,j,:)));
    else
    plot(squeeze(res.lfpresp(:,j,:)));
    end
    
end

function PlotByStim(res, rc, type)

[r,c] = Nsubplots(size(res.lfpresp,3));

crange = minmax(res.lfpresp(:,:,2:end));
for j = 1:size(res.lfpresp,3)
    subplot(r,c,j);
    if type == 1
    imagesc(squeeze(res.lfpresp(:,:,j)));
    if j > 1
    caxis(crange);
    end
    else
    plot(squeeze(res.lfpresp(:,:,j)));
    end
end

function res = PlotMovie(res, rc, delay, mtimes)

sustep = 1000;
lfpstep = 10000;
lfptimes = mtimes;
if res.lineplot == 3
    for j = 1:size(res.lfpresp,2)
      lfpresp(:,j,:) = CalcCSD(squeeze(res.lfpresp(:,j,res.probes))','smooth',res.csdsk)';
    end
    res.lfpresp = lfpresp;
    if isfield(res,'lfprespa')
        for j = 1:size(res.lfprespa,2)
            lfpresp(:,j,:) = CalcCSD(squeeze(res.lfprespa(:,j,res.probes))','smooth',res.csdsk)';
        end
        res.lfprespa = lfpresp;
    end
        
end

if isfield(res,'lfpresp')
    
res.lfprange(1) = min(res.lfpresp(:));
res.lfprange(2) = max(res.lfpresp(:));
end
res.mtimes = mtimes;
PlotTimeSlice(res, res.mtimes(1));


fb = findobj(gcf,'Tag','NextFrame');
if isempty(fb)
    bp = [10 10 40 20];
    uicontrol(gcf, 'style','pushbutton','string','+','Position',bp,...
        'Tag','NextFrame','Callback',{@NextFrame, 0});
end
fb = findobj(gcf,'Tag','LastFrame');
if isempty(fb)
    bp = [50 10 40 20];
    uicontrol(gcf, 'style','pushbutton','string','-','Position',bp,...
        'Tag','LastFrame','Callback',{@NextFrame, -1});
end
fb = findobj(gcf,'Tag','PlayFwd');
if isempty(fb)
    bp = [90 10 40 20];
    uicontrol(gcf, 'style','pushbutton','string','>>','Position',bp,...
        'Tag','PlayFwd','Callback',{@PlayFrames,  1});
end
fb = findobj(gcf,'Tag','PlayBack');
if isempty(fb)
    bp = [130 10 40 20];
    uicontrol(gcf, 'style','pushbutton','string','<<','Position',bp,...
        'Tag','PlayBack','Callback',{@PlayFrames, -1});
end
fb = findobj(gcf,'Tag','PlayStop');
if isempty(fb)
    bp = [170 10 60 20];
    uicontrol(gcf, 'style','check','string','Stop','Position',bp,...
        'Tag','PlayStop');
    bp(1) = bp(1)+bp(3);
    bp(3) = 200;
    res.slider = uicontrol(gcf, 'style','slider','string','t','Position',bp,...
        'Tag','PlaySlider','Min',0,'max',200,'value',0,'sliderstep',[0.005 0.025],'CallBack',@FrameSlider);
else
    set(fb(1),'value',0);
end
set(gcf,'UserData',res);
ofig = gcf;
subplot(1,1,1); hold off;

for ti = 1:length(res.mtimes)
    stop = get(findobj(gcf,'Tag','PlayStop'),'value');
    if stop
        ti = NaN;
    else


        set(0,'currentfigure',ofig);
    PlotTimeSlice(res,res.mtimes(ti));
    res.inow = ti;
    set(gcf,'UserData',res);
    drawnow;
    if delay > 1
        pause
    else
       pause(delay);
    end
    end
end

function PlotTimeSlice(res,t);
    
if isfield(res,'lfprespa')
    subplot(2,1,1);
    [a,b] = min(abs(t-res.lfptimes));
    if res.plot

    imagesc(res.x(:,1),res.probes,squeeze(res.lfprespa(b,:,:))' - squeeze(res.lfpresp(b,:,:))');
    else
        imagesc(res.x(:,1),res.probes,squeeze(res.lfprespa(b,:,:))');
    end
    caxis(res.lfprange);
    if isfield(res,'titlestr')
        title(sprintf('%.1f ms %s',res.lfptimes(b)./10,res.titlestr{2}));
    else
        title(sprintf('%.1f ms',res.lfptimes(b)./10));
    end
end

if isfield(res,'lfpresp')
    [a,b] = min(abs(t-res.lfptimes));
    if res.plot & isfield(res,'lfprespa')
        title(sprintf('%.1f ms',res.lfptimes(b)./10));
        imagesc(squeeze(res.lfprespa(b,:,:)+res.lfpresp(b,:,:))'./2);
    else
        imagesc(res.x(:,1),res.probes,squeeze(res.lfpresp(b,:,:))');
    end
    caxis(res.lfprange);
    if isfield(res,'titlestr')
        title(sprintf('%.1f ms %s',res.lfptimes(b)./10,res.titlesr{1}));
    else
        title(sprintf('%.1f ms',res.lfptimes(b)./10));
    end
    if isfield(res,'slider') && res.lfptimes(b) <= 2000
    set(res.slider,'value',res.lfptimes(b)./10);
    end
end

function FrameSlider(a,b,step)

D = get(gcf,'UserData');
t = get(a,'value') .* 10;
[dt, it] = min(abs(t-D.mtimes)); 
PlotTimeSlice(D,D.mtimes(it));
drawnow;
D.inow = it;
set(gcf,'UserData',D);

