function [resps, details] = AbsSlant(varargin)
% simulate slant responses of simple Gaussian rf with a disparity tuning
% fucntion.

iw = 256;
ih = iw;
im = zeros(iw,ih);
rfpos(1) = -2;
rfpos(2) = 0;
dtscale = [0 1]; %offset, baseline
addisps = [-0.2:0.1:0.2];
pix2deg = 0.05;  %20 pixels per degree
stimr = 3;
dtrange = [];

tags = {'DT', 'Pcolor', 'RF', 'TiltLines'};
dtype = 'gaboramp';

j = 1;
while j <= length(varargin)
    if isstruct(varargin{j})
        details = varargin{j};
        resps = varargin{j+1};
        j = j+1;
    elseif strncmpi(varargin{j},'aa',2)
    elseif strncmpi(varargin{j},'dtscale',3)
        j = j+1;
        dtscale = varargin{j};
    elseif strncmpi(varargin{j},'dxs',3)
        j = j+1;
        addisps = varargin{j};
    elseif strncmpi(varargin{j},'dtrange',3)
        j = j+1;
        dtrange = varargin{j};
    elseif strncmpi(varargin{j},'even',2)
        dtype = 'evengabor';
    elseif strncmpi(varargin{j},'plotall',3)
        GetFigure(tags{1});
        hold off;
        plot(details.dtx,details.dtf);
        hold on;
        plot(details.disps,max(details.dtf) .* ones(size(details.disps)),'v');
        hold off;
        GetFigure(tags{2});
        hold off;
        imagesc(resps);
                GetFigure(tags{3});
        hold off;
        imagesc(details.rf);
        hold on;
        a = 0:pi./20:2*pi;
        r = ones(size(a)) .* stimr./pix2deg;
        [x,y] =pol2cart(a,r);
        plot(x+iw/2,y+ih/2,'w--');
                    GetFigure(tags{4});
        hold off;
        plot(details.slants,resps);

    elseif strncmpi(varargin{j},'rfp',3)
        j = j+1;
        rfpos = varargin{j};
    end
    j = j+1;
end

[x,y] = meshgrid([1:iw],[1:ih]);
x = (x-iw/2) .* pix2deg;
y = (y-ih/2) .* pix2deg;
sigma = 1;
slant = 0;
dg=0;
details.x = unique(x(:));
details.y = unique(y(:));
details.iw = iw;
details.ih = ih;
details.stimr = stimr;
details.pix2deg = pix2deg;
dgs = 0:0.025:0.1;
slants = 0:pi/8:2*pi;
dg = 0.15;

r = (x-rfpos(1)).^2 + (y - rfpos(2)).^2;
rf = exp(-(r./(2.*sigma.^2)));
details.slants = slants;
details.disps = addisps;
[a,b] = min(abs(addisps - mean(addisps)));
midd = addisps(b);
%imagesc(rf);
for k = 1:length(slants)
    slant = slants(k);
    xo = x .* cos(slant) + y .* sin(slant);
    for j = 1:length(addisps)
        %    dg = dgs(j);
        disps = addisps(j) + xo .* dg;

        % now calculated
        pixresps = GetResps(disps,dtype,dtscale) .* rf;
        resps(j,k) = sum(pixresps(:))./sum(rf(:));
        if addisps(j) == midd;
            details.planes{k} = disps;
        end
    end
end
details.rf = rf;
if isempty(dtrange)
    details.dtx = [(min(details.planes{1}(:))-range(addisps)/2):0.01:(max(details.planes{end}(:))+range(addisps)/2)];
else
    details.dtx = dtrange;
end
details.dtf = GetResps(details.dtx,dtype,dtscale);



function resps = GetResps(disps, type, dtscale)

if strncmpi(type,'linear',3)
    resps = DispResp(disps, type);
elseif strncmpi(type,'gaboramp',8)
    resps = DispResp(disps,type,[0.01 0.3 -pi/2 10 0 0],0.8);
elseif strncmpi(type,'evengabor',3)
    resps = DispResp(disps,'gabor',[0.01 0.2 0 10 0.2 0],1);
elseif strncmpi(type,'gabor',3)
    resps = DispResp(disps,type,[0.01 0.2 -pi/2 10 0 0],1);
elseif strncmpi(type,'ramp',3)
    resps = DispResp(disps, type);
else
    resps = DispResp(disps, type);
end

if nargin > 2
    resps = (resps+dtscale(1)) .* dtscale(2);
end

function resps = DispResp(disps,type, varargin)
%apply disparity tuning curve to a set of disparities.


if strncmpi(type,'linear',3)
resps = disps .* 10;
elseif strncmpi(type,'gabor',3)
    sd = varargin{1}(2);
    k = varargin{2};
    resps = k .* Gabor(varargin{1},'xv',disps,'pix2deg',0.05) + (1-k) .* erf(disps./sd);
elseif strncmpi(type,'ramp',3)
    resps = disps.*10;
    resps(find(resps>10)) = 10;
    resps(find(resps < -10)) = -10;
else
    resps = disps .* 10;
end