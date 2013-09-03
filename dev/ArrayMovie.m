function ArrayMovie(Expt, ids, varargin)
%Plot an array Expt data strcture as a movie

timerange = [];
j = 1; 
while j <= length(varargin)
    if strncmpi(varargin{j},'timerange',6)
        j = j+1;
        timerange = varargin{j};
    end
    j = j+1;
end
if isfield(Expt,'Header')
    tid = find(ismember([Expt.Trials.id],ids));
    res.lfpresp = mean(cat(3,Expt.Trials(tid).LFP),3);
    Array = GetArrayConfig(Expt.Header.loadname);
    res.lfptimes = Expt.Header.lfptimes;
elseif isnumeric(Expt)
    res.lfpresp = Expt;
    if isfield(ids,'sdfs') && isfield(ids.sdfs,'lfptimes') %RevCor return struct
        Array = GetArrayConfig(ids.Header.loadname);
        res.lfptimes = ids.sdfs.lfptimes;
    else
        Array = GetArrayConfig(ids.loadname);
        res.lfptimes = ids.times;
    end
end

if ~isempty(timerange)
    id = find(res.lfptimes > timerange(1) & res.lfptimes < timerange(end));
    res.lfptimes = res.lfptimes(id);
    res.lfpresp = res.lfpresp(id,:);
end
    
nx = max(Array.X);
ny = max(Array.Y);
res.arraysize = [nx ny];
res.probemap = sub2ind([nx ny],Array.X,Array.Y);
PlotMovie(res,0.01, res.lfptimes);

function NextFrame(a,b,step)

D = get(gcf,'UserData');
t = D.inow;
if step == 0
    t = t+1;
    if t > length(D.mtimes)
        D.mtimes = [D.mtimes D.mtimes(end)+mean(diff(D.mtimes))];
        t = length(D.mtimes);
    end
elseif step < 0
    t = t+step;
    if t < 1
        t = 1;
    end
end
PlotTimeSlice(D,D.mtimes(t));   
D.inow = t;
set(gcf,'UserData',D);

function FrameSlider(a,b,step)

D = get(gcf,'UserData');
t = get(a,'value');
it = round(t .* length(D.mtimes)./200);
if it > 0
    PlotTimeSlice(D,D.mtimes(it));
end
drawnow;
D.inow = it;
set(gcf,'UserData',D);

function PlayFrames(a,b,step)

D = get(gcf,'UserData');
t = D.inow;
set(findobj(gcf,'Tag','PlayStop'),'value',0);
if step > 0
    if D.inow >= length(D.mtimes)
        D.inow = 1;
    end
    t = [D.inow:step:length(D.mtimes)];
elseif step < 0
    t = [D.inow:step:1];
end
for it = t
    stop = get(findobj(gcf,'Tag','PlayStop'),'value');
    if stop
        break;
    end
    PlotTimeSlice(D,D.mtimes(it));
    drawnow;
end
D.inow = it;
set(gcf,'UserData',D);

function res = PlotMovie(res, delay, mtimes)

sustep = 1000;
lfpstep = 10000;
    


if isfield(res,'spkresps')
res.surange(1) = min(res.spkresps(:));
res.surange(2) = max(res.spkresps(:));
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
subplot(2,1,1); hold off;
subplot(2,1,2); hold off;
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

nplots = 1;
if isfield(res,'lfprespa')
    subplot(2,1,1);
    nplots = 2;
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
elseif isfield(res,'spkresps')
    subplot(2,1,1);
    nplots = 2;
    [a,b] = min(abs(t-mean(res.times,1)));
    imagesc(squeeze(res.spkresps(b,:,:))');
    caxis(res.surange);
    title(sprintf('%.1f ms',res.times(1,b)./10));
end

if isfield(res,'lfpresp')
    if nplots == 1
        subplot(1,1,1);
    else
        subplot(2,1,2);
    end
    [a,b] = min(abs(t-res.lfptimes));
    if isfield(res,'lfprespa') & res.plot
        title(sprintf('%.1f ms',res.lfptimes(b)./10));
        imagesc(squeeze(res.lfprespa(b,:,:)+res.lfpresp(b,:,:))'./2);
    elseif isfield(res,'probemap')
        im = zeros(res.arraysize);
        im(res.probemap) = res.lfpresp(b,:);
        imagesc(im);
    else
        imagesc(res.x(:,1),res.probes,squeeze(res.lfpresp(b,:,:))');
    end
    caxis(res.lfprange);
    if isfield(res,'titlestr')
        title(sprintf('%.1f ms %s',res.lfptimes(b)./10,res.titlestr{1}));
    else
        title(sprintf('%.1f ms',res.lfptimes(b)./10));
    end
    if isfield(res,'slider')
        t = res.lfptimes(b)./10;
        if t <= get(res.slider,'Max') & t >= get(res.slider,'Min')
            set(res.slider,'value',res.lfptimes(b)./10);
        end
    end
end


function PlotTimeCourseProbes(rc, probes, style, sumy)

showvar = 1;

[nr,nc] = Nsubplots(length(probes));
if sumy == 2
    ny = ceil(length(rc.yvals)/2);
    yv = rc.yvals(1:ny);
    zv = rc.yvals((ny+2):length(rc.yvals));
end

cr = [min(rc.lfp(:)) max(rc.lfp(:))];
for probe = probes
    subplot(nr,nc,probe);
    if sumy == 1
        if style == 1
            plot(rc.lfptimes./10,rc.lfpresp(:,:,probe)');
            if showvar
                v = var(rc.lfpresp(:,:,probe),[],2);
                scale = max(max(rc.lfpresp(:,:,probe)))./max(v);
                hold on;
                plot(rc.lfptimes./10,v.*scale,'k');
                hold off;
                end
        else
            imagesc(rc.x(:,1),rc.lfptimes./10,rc.lfpresp(:,:,probe));
            caxis(cr);
        end
    elseif sumy == 2  %diff
        if style == 1
            plot(rc.lfptimes./10,sum(rc.lfp(:,:,yv,probe),3)-sum(rc.lfp(:,:,zv,probe),3));
        else
            imagesc(rc.x(:,1),rc.lfptimes./10,sum(rc.lfp(:,:,yv,probe),3)-sum(rc.lfp(:,:,zv,probe),3));
            caxis(cr);
        end
    end
end
