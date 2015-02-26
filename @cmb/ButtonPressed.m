function ButtonPressed(src, data)
global mousept;

DATA = GetDataFromFig(src);

if isempty(DATA) %broken - lost parent figure?
    return;
end
it = findobj(src,'Tag','Clusterid');
if ~isempty(it)
    mousept.cluster = get(it,'value');
else
    mousept.cluster = 1;
end

start = get(gca,'CurrentPoint');
xl = get(gca,'xlim');
yl = get(gca,'ylim');
if start(1,1) < xl(1) || start(1,1) > xl(2) || start(1,2) < yl(1) || start(1,2) > yl(2)
    return; %out of range
end
mousept.drags = 0;

mousept.color = DATA.spkcolor{mousept.cluster+1};
mousept.lasth = 0;
lastkey = get(gcf,'CurrentCharacter');

% don't want the axis rescaling as we draw the ellipse
set(gca,'Xlimmode','Manual','Ylimmode','Manual');
mousept.down = 1;
if isfield(mousept,'mode')
    oldmode = mousept.mode;
else
    oldmode = 0;
end
mousept.mode = strmatch(get(gcf,'SelectionType'),{'normal' 'alternate'  'extend'  'open'});
if oldmode == 10  && mousept.mode ==1% a left button press after scrolling
    mousept.mode = 10;
    return;
end
p = get(gca,'UserData');
mousept.xrange = diff(get(gca,'Xlim'));
mousept.yrange = diff(get(gca,'Ylim'));

if isempty(p)
    p = DATA.probe;
end
if mousept.cluster <= size(DATA.cluster,1) && p <= size(DATA.cluster,2)
    C = DATA.cluster{mousept.cluster,p};
else
    C.h = 0;
end
mousept;
hold on; %othewise drawing ellipse deletes data
if mousept.mode == 1
    mousept.start = get(gca,'CurrentPoint');
    
    if size(DATA.cluster,1) >= mousept.cluster & size(DATA.cluster,2) >= p &...
            isfield(DATA.cluster{mousept.cluster,p},'h')
        for j = 1:size(DATA.cluster,1)
            distance(j) = cmb.DistanceToCluster(DATA.cluster{j,p},mousept.start(1,1:2));
        end
        if(min(distance) > 1.05) % pressed outside = start over
            if ishandle(C.h)
                delete(C.h);
            end
        else  %% pressed inside; select this cluster, move if mouse moves
            [d, cl]= min(distance);
            C = DATA.cluster{cl,p};
            if ~isempty(it)
                set(it,'value',cl);
            end
            fi = findobj(src,'Tag','Clusterid');
            if isfield(C,'forceid') && C.forceid > 0
                DATA.forceclusterid = C.forceid;
                if ~isempty(fi)
                    set(fi,'value',C.forceid);
                end
            else
                DATA.forceclusterid = 0;
                if ~isempty(fi)
                    set(fi,'value',cl);
                end
            end
            DATA.cluster{cl,p}.h = cmb.DrawCluster(DATA.cluster{cl,p}, DATA.spkcolor{cl+1});
            mousept.mode = 5;
            mousept.cluster = cl;
            mousept.color = DATA.spkcolor{mousept.cluster+1};
            mousept.c = [DATA.cluster{cl,p}.x(1) DATA.cluster{cl,p}.y(1)];
            mousept.r = [DATA.cluster{cl,p}.x(2) DATA.cluster{cl,p}.y(2)];
            DATA.cluster{cl,p}.touched = 1;
            mousept.offset = mousept.start(1,1:2) - mousept.c;
            mousept.angle = -DATA.cluster{cl,p}.angle;
            mousept.lasth = DATA.cluster{cl,p}.h;
            set(DATA.cluster{cl,p}.h,'linewidth',2);
            if oldmode == 5 %second press in ellipse - move it
                mousept.down = 1;
                mousept.mode = 11;
                mousept.lasth = DATA.cluster{cl,p}.h;
            else
                mousept.down = 1;  %%ignore drag, release
            end
            set(DATA.toplevel,'UserData',DATA);
        end
    else
        p
    end
elseif mousept.mode == 2 %R button
    if isfield(C,'h') & C.h & ishandle(C.h)
        delete(C.h);
    end
    cl = mousept.cluster;
    mousept.c = [DATA.cluster{cl,p}.x(1) DATA.cluster{cl,p}.y(1)];
    mousept.r = [DATA.cluster{cl,p}.x(2) DATA.cluster{cl,p}.y(2)];
    mousept.angle = -DATA.cluster{cl,p}.angle;
    mousept
    mousept.start = get(gca,'CurrentPoint');
    if mousept.start(1,1) > mousept.c(1) + mousept.r(1)
        mousept.mode = 6;
    elseif mousept.start(1,1) < mousept.c(1) - mousept.r(1)
        mousept.mode = 7;
    elseif mousept.start(2,2) > mousept.c(2) + mousept.r(2)
        mousept.mode = 8;
    elseif mousept.start(2,2) < mousept.c(2) - mousept.r(2)
        mousept.mode = 9;
    end
    fprintf('Start %.2f,%.2f C %.2f,%.2f, R%.2f %.2f, mode %d\n',...
        mousept.start(1,1),mousept.start(2,2),mousept.c(1),mousept.c(2),mousept.r(1),mousept.r(2),mousept.mode);
elseif mousept.mode == 3 % R button
    if ishandle(C.h) delete(C.h); end
    
end


