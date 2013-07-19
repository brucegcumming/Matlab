function [allresps, details] = dori(varargin)
%simulate a 2-D BEM and test with slanted stimuli etc.
%
%
ori = 0;
dori = 0;
dsf = 0;
sf = 1;
sd = [0.4 0.39];
nruns = 1000;
disp = 4;
w = 256;
h = 256;
ndots = 10000;
disps = [-0.3:0.05:0.3];
dgs = [-0.6:0.1:0.6];
dfs = [0 0 0];
dfs = zeros(size(dgs));
calcbw = 0;
silent = 0;
j = 1;
while j <= length(varargin) 
    if strncmpi(varargin{j},'bw',2)
        calcbw = 1;
    elseif strncmpi(varargin{j},'nruns',3)
        j = j+1;
        nruns = varargin{j};
    elseif strncmpi(varargin{j},'ori',3)
        j = j+1;
        ori = varargin{j};
    elseif strncmpi(varargin{j},'dori',3) %RF dori
        j = j+1;
        dori = varargin{j};
    elseif strncmpi(varargin{j},'disps',3) %RF dori
        j = j+1;
        disps = varargin{j};
    elseif strncmpi(varargin{j},'dsf',3) %RF dori
        j = j+1;
        dsf = varargin{j};
    elseif strncmpi(varargin{j},'dgs',3)
        dgs = [-0.3 -0.2 -0.1 0 0.1 0.2 0.3];
        if length(varargin)  > j & isnumeric(varargin{j+1})
            j = j+1;
            dgs = varargin{j};
        end
        dfs = zeros(size(dgs));
    elseif strncmpi(varargin{j},'dfs',3)
        dfs = [-0.3 -0.2 -0.1 0 0.1 0.2 0.3];
        if length(varargin)  > j & isnumeric(varargin{j+1})
            j = j+1;
            dfs = varargin{j};
        end
        dgs = zeros(size(dfs));
    elseif strncmpi(varargin{j},'ndots',3) %RF dori
        j = j+1;
        ndots = varargin{j};
    elseif strncmpi(varargin{j},'silent',3) %RF dori
        silent =1;
    elseif strncmpi(varargin{j},'track',3) %RF dori
        silent =2;
        
    elseif strncmpi(varargin{j},'seed',3) %RF dori
        j = j+1;
        seed = varargin{j};
        rand('twister',seed);
    end
    j = j+1;
end

details.dgs = dgs;
details.dfs = dfs;
degpix = 0.02;

lrfa = Gabor([sf-dsf/2 sd(1) 0 1 0 0 ori-dori/2 sd(2)],'pix',degpix);
lrfb = Gabor([sf-dsf/2 sd(1) pi/2 1 0 0 ori-dori/2 sd(2)],'pix',degpix);
rrfa = Gabor([sf+dsf/2 sd(1) 0 1 0 0 ori+dori/2 sd(2)],'pix',degpix);
rrfb = Gabor([sf+dsf/2 sd(1) pi/2 1 0 0 ori+dori/2 sd(2)],'pix',degpix);
dmax = max(abs(disps));
colors = mycolors;

if calcbw
    sfs = 0.1:0.1:10;
    oris = ori-pi/2:pi/20:ori+pi/2;
    [X,Y] = meshgrid([1:w]-w/2,[1:h]-h/2);
for j = 1:length(sfs)
    sf = sfs(j);
    im = sin((X.*cos(ori) + Y .* sin(ori)) .* 2 .* pi .*sf .* degpix);
    allresps(j,1) = sum(lrfa(:) .* im(:)).^2 + sum(lrfb(:) .* im(:)).^2;
end
    [maxr, opt] = max(allresps);
    sf(1) = interp1(allresps(1:opt),sfs(1:opt),maxr/2);
    sf(2) = interp1(allresps(opt:end),sfs(opt:end),maxr/2)
    sfopt = sfs(opt);
    details.sfbw = sf(2)./sf(1);
for j = 1:length(oris)
    sf = sfs(j);
    im = sin((X.*cos(oris(j)) + Y .* sin(oris(j))) .* 2 .* pi .*sfopt .* degpix);
    oresps(j) = sum(lrfa(:) .* im(:)).^2 + sum(lrfb(:) .* im(:)).^2;
end
    [maxr, opt] = max(oresps);
    or(1) = interp1(oresps(1:opt),oris(1:opt),maxr/2);
    or(2) = interp1(oresps(opt:end),oris(opt:end),maxr/2);
    or = or .* 180./pi
    details.orbw = or(2) - or(1);
    return;
end

if silent == 0

GetFigure('Doria');
hold off;
end
[a, findk] = min(abs(disps - 0.2));
for m = 1:length(dgs)
%    divide dgs by 2, because add to each eye
    dg = dgs(m)/2;
    df = dfs(m)/2;
for k = 1:length(disps)
    disp = disps(k) ./ (2 .* degpix);
    if silent > 1
        fprintf('dx %.3f dg %.3f/%.3f at %s\n',disps(k),dg,df,datestr(now));
    end
    for j = nruns:-1:1
        iml = zeros(w,h);
        imr = iml;
        if ndots == 0
            spc = 0.025:0.05:1;
            x = [spc  ones(size(spc)) .*0.5];
            y = [ones(size(spc)) .* 0.5 spc];
            x = ceil(x.*w);
            y = ceil(y.*h);
        else
        x = ceil(rand(1,ndots) .* (w));
        y = ceil(rand(1,ndots) .* (h));
        end
        xr = round(x+disp+dg.*(y-h/2)+df.*(x-w/2));
        xl = round(x-disp-dg.*(y-h/2)-df.*(x-w/2));
        rid = find(xr > 0 & xr <= w);
        lid = find(xl > 0 & xl <= w);
        lix  = sub2ind(size(iml),xl(lid),y(lid));
        iml(lix) = 1;
        rix  = sub2ind(size(iml),xr(rid),y(rid));
        imr(rix) = 1;
        imr = imr - mean(imr(:));
        iml = iml - mean(iml(:));
        ba(j) = sum(iml(:).*lrfa(:)) .* sum(imr(:) .* rrfa(:));
        bb(j) = sum(iml(:).*lrfb(:)) .* sum(imr(:) .* rrfb(:));
    end
    if k == findk
        rim = rrfa;
        rim(rix) = -1;
        lim = lrfa;
        lim(lix) = -1;
    end
        alldisps(m,k) = disp;
        alldgs(m,k) = dg;
        alldfs(m,k) = df;
    resp(k) = mean(ba+bb);
    noises(k) = std(ba+bb);
end
if silent == 0 
    plot(disps,resp,'color',colors{m});
hold on;
end
allresps(:,m) = resp;
allsd(:,m) = noises;

end

details.allsd = allsd;
if silent == 0 

GetFigure('Dorib');
imagesc(allresps);
if length(dgs) ==1
GetFigure('Doric');
subplot(1,2,1); imagesc(lim')
subplot(1,2,2); imagesc(rim')
end
end
