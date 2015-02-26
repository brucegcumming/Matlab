function MDS = DavidShape(C, varargin)

v2d = 7;  %rough guees for converting unit distace in MDS to 50microns

v2mu = v2d * 50;


for p = 1:length(C)
    ms = C{p}.MeanSpike.ms;
    V = ms(p,:);
    amp =  V * ms';
    amp = amp ./amp(p);
    amps(p,:) = amp;
    sigma = std(ms,[],2);
    [a,b] = sort(sigma);
    stds(p,:) = sigma./sigma(p);
    MDS.mahal(p,:) = C{p}.mahal;
    MDS.meanspike{p} = C{p}.MeanSpike.ms;
    if ~isempty(C{p}.next)
        nc = 1;
        for j = 1:length(C{p}.next)
            nc = nc+1;
            ms = C{p}.MeanSpike.ms;
            MDS.xspike{p}.next{j}.meanspike = ms;
            V = ms(p,:);
            amp =  V * ms';
            amp = amp ./amp(p);
            amps(p,:) = amps(p,:) + amp;
            sigma = std(ms,[],2);
            [a,b] = sort(sigma);
            stds(p,:) = stds(p,:) + sigma'./sigma(p);
        end
        amps(p,:) = amps(p,:) /nc;
        stds(p,:) = stds(p,:) /nc;
    end
end

GetFigure('Amp');
imagesc(amps);
MDS.amps = amps;
MDS.sds = stds;
for p = 1:length(C)
for k = 1:p-1
    x = mean([amps(p,k) amps(k,p)]);
    if x > 1
        x = 0.95;
    end
    amps(p,k) = x;
    amps(k,p) = x;
end
end
mds = mdscale(amps,2);
GetFigure('MDS');
MDS.v2mu = v2mu;
gid = find(MDS.mahal(:,4) > 2);
bid = find(MDS.mahal(:,4) <= 2);
hold off;
plot(mds(gid,1).*v2mu,mds(gid,2).*v2mu,'ro','Buttondown',@HitMDS);
hold on;
plot(mds(bid,1).*v2mu,mds(bid,2).*v2mu,'o','Buttondown',@HitMDS);

for j = 1:length(mds)
    text(mds(j,1).*v2mu,mds(j,2).*v2mu,sprintf('%d',j),...
        'horizontalalignment','left','verticalalignment','bottom');
end
MDS.x = mds(:,1);
MDS.y = mds(:,2);
setappdata(gcf,'MDS',MDS);
X = mds(:,1)+i*mds(:,2);
d = [];
for p = 1:length(C)
    for k = 1:p-1
        d(end+1) = abs(X(p)-X(k));
        a(length(d)) = amps(p,k);
        ps(length(d),:) = [p k];
    end
    ds = abs(abs(X(p)-X));
    [x,b] = sort(MDS.amps(p,:),'descend');
    MDS.chspk(p,:) = b(1:4);
    [x,b] = sort(ds,'ascend');
    MDS.dspk(p,:) = b(1:4);
end
DATA.d = d;
DATA.amp = a;
DATA.probes = ps;
GetFigure('Distance-amp');
hold off;


setappdata(gcf,'Probes',DATA);

plot(d .* v2mu,a,'o','Buttondown',@HitScatter);
hold on;


%
%Exponential fit is approxiamte. Even a true exponential decay 
%gets distorted by discretized sampling, since the true peak may fall 
%between two probes. see matlab/sims/DistanceVoltage
dx = 0:0.01:1;
A = exp(-dx * v2d);
plot(dx .* v2mu,A,'k-','linewidth',2);
A = exp(-dx *(v2d * 2));
plot(dx .* v2mu,A,'r-','linewidth',1);
A = exp(-dx * (v2d / 2));
plot(dx .* v2mu,A,'r-','linewidth',1);
legend({'' 'Fit from laminar data' '*2 /2'});
xlabel('Microns');

function HitMDS(a,b,varargin)

x = get(gca,'currentpoint');
MDS = getappdata(gcf,'MDS');
x = x./MDS.v2mu;
[a,p] = min(abs(x(1,1)-MDS.x+i*(x(1,2)-MDS.y)));
fprintf('Probe %d\n',p);
GetFigure('Size');
hold off;
plot(MDS.amps(p,:));
hold on;
plot(MDS.sds(p,:),'r');
[a,b] = sort(MDS.amps(p,:),'descend');
for j = 1:6
    text(b(j),a(j),sprintf('%d',b(j)));
end
GetFigure('MeanSpike');
imagesc(MDS.meanspike{p});




function HitScatter(a,b,varargin)

p = get(gca,'currentpoint');
P = getappdata(gcf,'Probes');
[a,b] = min(abs(p(1,1)-P.d+i*(p(1,2)-P.amp)));
fprintf('Probes%s\n',sprintf(' %d',P.probes(b,:)));

