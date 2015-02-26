function [x,y,z] = CalcDensity(DATA, expspks, mode)
%[x,y,z] = CalcDensity(DATA, expspks, mode)
%Calcuate density plot for spike properties in Combine

if ismember(DATA.syncsign,[-1 1])
energy = DATA.Spikes.cx(DATA.sids{1});
vw = DATA.Spikes.cy(DATA.sids{1});
elseif isfield(DATA,'AllClusters')
    if iscell(DATA.AllClusters)
    energy  = DATA.AllClusters{DATA.currentexpt}(DATA.probe).cx(expspks);
    vw  = DATA.AllClusters{DATA.currentexpt}(DATA.probe).cy(expspks);
    else
    energy  = DATA.AllClusters(DATA.probe).cx(expspks);
    vw  = DATA.AllClusters(DATA.probe).cy(expspks);
    end
else
energy = DATA.Spikes.cx(expspks);
vw = DATA.Spikes.cy(expspks);
end
DATA.plot.DensitySigma = [3 3];

if length(vw) > 10000
    lprc = 0.01;
    hprc = 99.99;
elseif length(vw) > 10000
    lprc = 0.1;
    hprc = 99.9;
    
elseif length(vw) > 1000
    lprc = 1;
    hprc = 99;
    DATA.plot.DensitySigma = [5 5];
else
    lprc = 5;
    hprc = 95;
    DATA.plot.DensitySigma = [10 10];
end    

if DATA.xyfig == gcf && mode < 3
erange = get(gca,'Xlim');
vrange = get(gca,'Ylim');
else
    if prctile(energy,hprc) > prctile(energy,hprc-1) * 3
        erange = [prctile(energy,lprc) prctile(energy,hprc-1)];
    else
        erange = [prctile(energy,lprc) prctile(energy,hprc)];
    end
vrange = [prctile(vw,lprc) prctile(vw,hprc)];
end
%GetFigure('DensityPlot');
hold off;
nbins = 200;
[x,y] = meshgrid(linspace(erange(1),erange(2),nbins),linspace(vrange(1),vrange(2),nbins));
tic;
if mode ==1 % add real gaussian to grid for each
    z = zeros(size(x));
    sx = (diff(erange)/100)^2;
    sy = (diff(vrange)/100)^2;
    for j=1:length(energy)
        z = z + exp(-(x-energy(j)).^2/sx - (y-vw(j)).^2/sy);
    end
elseif mode ==2 || mode == 3 %build fine 2-D histogram, then smooth
    [gx,gy] = meshgrid(-10:10,-10:10);
    sx= DATA.plot.DensitySigma(1);
    sy= DATA.plot.DensitySigma(2);
    G = exp(-(gx).^2/sx - (gy).^2/sy);
    G = G./sum(G(:));
    z = zeros(size(x));
% ignore spikes where both are set to 0
    idx = find(vw ~=0 | energy ~=0);
    
    vi = 1+floor(nbins * (vw(idx)-vrange(1))/diff(vrange));
    ei = 1+floor(nbins * (energy(idx)-erange(1))/diff(erange));
    idx = find(ei > 0 & ei <= nbins & vi > 0 & vi <= nbins);
    for j =idx
        z(vi(j),ei(j)) = z(vi(j),ei(j))+1;
    end
    z = conv2(z,G,'same');
end
