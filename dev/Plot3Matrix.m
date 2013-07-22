function Plot3Matrix(X, varargin)

res.X = X;
mtimes = 1:size(X,1);
res.range(1) = min(X(:));
res.range(2) = max(X(:));
res.mtimes = mtimes;
PlotMovie(res,0, mtimes);

function res = PlotMovie(res, delay, mtimes)

    
PlotTimeSlice(res, res.mtimes(1),1);


fb = findobj(gcf,'Tag','NextFrame');
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
    PlotTimeSlice(res,res.mtimes(ti),1);
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

function PlotTimeSlice(res,t,doset);
    

[a,b] = min(abs(t-res.mtimes));
    imagesc(squeeze(res.X(b,:,:))');
    caxis(res.range);
    if isfield(res,'titlestr')
        title(sprintf('%.1f ms %s',res.mtimes(b)./10,res.titlesr{1}));
    else
        title(sprintf('%.1f ms',res.mtimes(b)./10));
    end
    if doset & isfield(res,'slider') && res.mtimes(b) <= 2000
        t = t .* 200./res.mtimes(end);
    set(res.slider,'value',t);
    end

function FrameSlider(a,b,step)

D = get(gcf,'UserData');
t = get(a,'value');
[dt, it] = min(abs((t .* D.mtimes(end)./200)-D.mtimes)); 
PlotTimeSlice(D,D.mtimes(it),0);
drawnow;
D.inow = it;
set(gcf,'UserData',D);

