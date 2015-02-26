function WorldSeries(ps, varargin)

if nargin == 0
    ps = [0.5:0.02:0.64  0.65:0.005:0.7 0.72:0.02:1];
end
nrpt = 1000;
nresults = 20; 

s = scanlines('WorldSeriesData.txt');
ng = 0;
for j = 1:length(s)
    id = regexp(s{j},'[a-z]\s+4');
    if ~isempty(id)
        start = id(1)+1;
        while s{j}(start) ~= '4'
            start = start+1;
        end
        start = start+2;
        ng = ng+1;
        [x,y] = sscanf(s{j}(start:end),'%f');
        games(ng) = 4+x;
    end
end
truemean = mean(games);
nresults = ng;

gsum=0;
for j = length(games):-1:1
    gsum = gsum+games(j);
end

for j = 1:length(ps)
    p = ps(j);
    for l = 1:nrpt
    wins = binornd(1,p,7,nresults);
    wa = cumsum(wins ==1);
    wb = cumsum(wins ==0);
    for k = 1:size(wa,2)
        a =find(wa(:,k) ==4 | wb(:,k) ==4);
        a = min(a);
        [lens(k,l)] = a;
        winner(k,l) = wa(a,k) ==4;
    end
    end
    m = mean(lens);
    w = mean(winner);
    slen(j) = mean(m); %mean length
    sds(j) = std(m);
    winrate(j) = mean(w);
end
if length(ps) > 1
    errorbar(winrate, slen,2.*sds,'o-');
    hold on;
    plot(minmax(ps),[truemean truemean],'r-');    
else
    hist(lens(:));
end