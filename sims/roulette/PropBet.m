function [pots, all] = PropBet(edge, nruns, varargin)

%
% need to keep track of nwins/ /nbets to debug zero values.

payout = 6;
pwin = edge/6;
runl = 1000;
bets = [1 2];
all.wins = [];
maxl = 1000;
seed = 1;


j = 1;
while j < nargin-1
    if strncmpi(varargin{j},'len',3)
        runl = varargin{j+1};
        j = j+1;
    elseif strncmpi(varargin{j},'bets',3)
        bets = varargin{j+1};
        j = j+1;
    elseif strncmpi(varargin{j},'maxl',4)
        maxl = varargin{j+1};
        j = j+1;
    elseif strncmpi(varargin{j},'seed',4)
        seed = varargin{j+1};
        j = j+1;
    end
    j = j+1;
end
rand('state',seed);


%all = CalcWinsOnOff(nruns, runl, 1./6, 6, 1);
pwin = (1 + edge)./6;
[all, pots] = CalcWinsStop(nruns, runl, pwin, 6, 1, maxl);
%CalcWinsProp(nruns, runl, bets);

function [all, totals] = CalcWinsStop(nruns, runl, pwin, payout, bets, maxl)

figa = 'OneRun';
figb = 'Summary';
%spins = rand(nruns,runl);
GetFigure(figa);
subplot(1,1,1);
savedata = 0;
plotrun = 0;

for j = nruns:-1:1;
    spins = rand(1,runl);
    wins = spins < pwin;
    pot = 0.001;
    tpot = pot;
    bet = 1;
    lastbet = bet;
    nwin = [0 0];
    ncross = 0;
    lastpot = pot;
    didbet = zeros(1,runl);
    pots(runl) = 0;
    tpots(runl) = 0;
    all.wins(runl) = 0;
    allbets(runl) = 0;
    if savedata
    alltpots(j,runl) = 0;
    end
    for k = 1:runl;
        if tpot < 0
            didbet(k) = 0;
            if wins(k);
                tpot = tpot + bet * payout;
            end
        else
            didbet(k) = 1;
            pot = pot - bet;
            if wins(k);
                tpot = tpot + bet * payout;
                pot = pot + bet * payout;
            end
        end
        tpot = tpot - bet;
        all.wins(k) = sum(nwin);
        pots(k) = pot;
        tpots(k) = tpot;
        allbets(k) = bet;
        if savedata
        alltpots(j,k) = tpot;
        end
    end
    sgn = diff(sign(tpots))/2;
    id = find(abs(sgn)>0);
    if tpots(end) < 0 %finish negative 
        sgn = [sgn 1];
    else
        sgn = [sgn -1];
    end
%inlcude end
    id = [id runl];
    xids{j} = id.*sgn(id); %list of crossings
    if plotrun > 1
    hold off; plot(tpots);
    end
        if j == 1
        bets = cumsum(didbet);
        end
        bets = cumsum(didbet);
    if length(id) > 1
        for m = 1:length(id)
            nbets{j}(m) = bets(id(m));
        end
        lid = find(diff(id) > maxl & sgn(id(2:end)) > 0);
        if length(lid)
            stops(j) = id(lid(1))+1; %% first netagive tpot of run
            stoppots(j) = pots(stops(j));
            allstoppots{j} = pots(id);
            if plotrun
            hold off; plot(tpots);
            hold on;
            plot([id(lid(1)) id(lid(1)+1)],[0 0],'-','linewidth',3);
            title(sprintf('p %.3f run %d',pwin,j));
            drawnow;
            end
        else
            stoppots(j) = pots(id(m));
            allstoppots{j} = pots(id);
            stops(j) = NaN; %% first netagive tpot of run
        end
    else
            stoppots(j) = pots(end);
            stops(j) = NaN;
            allstoppots{j} = pots(end);
        for m = 1:length(id)
            nbets{j}(m) = bets(id(m));
        end
    end
    ttotals(j) = tpots(k);
    totals(j) = pots(k);
end
GetFigure(figb);
subplot(2,2,1);
hold off;
plot(stops,stoppots,'o');

all.stops = stops;
all.stoppots = stoppots;
all.xids = xids;
all.nbets = nbets;
all.allstoppots = allstoppots;
if savedata
all.data = alltpots;
end
function all = CalcWinsOnOff(nruns, runl, pwin, payout, bets)

spins = rand(nruns,runl);
wins = spins < pwin;
for j = nruns:-1:1;
    pot = 0.001;
    tpot = pot;
    bet = 1;
    lastbet = bet;
    nwin = [0 0];
    ncross = 0;
    lastpot = pot;
    for k = runl:-1:1;
        tpot = tpot - bet;
        if tpot < 0
            didbet(j,k) = 0;
            if wins(j,k);
                tpot = tpot + bet * payout;
            end
        else
            didbet(j,k) = 1;
            pot = pot - bet;
            if wins(j,k);
                tpot = tpot + bet * payout;
                pot = pot + bet * payout;
            end
        end
        all.wins(j,k) = sum(nwin);
        pots(j,k) = pot;
        tpots(j,k) = tpot;
        allbets(j,k) = bet;
    end
end
subplot(2,2,1);
hold off;
plot(tpots(:,1),pots(:,1),'o');
refline(1);
hold on;
[h,p] = ttest(pots(:,1));
fprintf('Mean Winnigs %.3f +- %.3f(p < %.6f)\n',mean(pots(:,1)),std(pots(:,1)),p);
fprintf('Theoretical Winnings %.3f +- %.3f\n',mean(tpots(:,1)),std(tpots(:,1)));
ylabel('acutal gain');
xlabel('theoretical gain');
subplot(2,2,3);
hist(tpots(:,1),25);
subplot(2,2,2);
hist(pots(:,1),25);

all.payout = payout;
all.betvals = bets;
all.bets = allbets;
all.pots = pots;
all.tpots = tpots;
all.nbets = sum(didbet,2);

function CalcWinsProp(nruns, runl, bets)

spins = rand(nruns,runl);
for j = 1:nruns;
    pot = 0.001;
    bet = bets(1);
    lastbet = bet;
    nwin = [0 0];
    ncross = 0;
    lastpot = pot;
    for k = 1:runl;
        if pot > 0
            bet = bets(2);
            bid = 2;
        else
            bet = bets(1);
            bid = 1;
        end
        if sign(lastpot) ~= sign(pot)
            ncross = ncross + 1;
        end
        lastpot = pot;
        if spins(j,k) < pwin %win
            pot = pot + bet * (payout-1);
            nwin(bid) = nwin(bid) + 1;
        else
            pot = pot - bet;
        end
        all.wins(j,k) = sum(nwin);
        pots(j,k) = pot;
        allbets(j,k) = bet;
    end
    all.ppot(j) = nwin(1) * bets(1) * payout + nwin(2) * bets(2) * payout - sum(allbets(j,:) < bets(2)) * bets(1)  - sum(allbets(j,:) > bets(1)) * bets(2);
    meanbets(j) = mean(allbets(j,:));
    crosses(j) = ncross;
    evenpay(j) = (sum(nwin) * (payout-1) - (runl-sum(nwin))) * meanbets(j);
    nwins(j,:) = nwin;
end
subplot(2,2,1);
hold off;
plot(pots(:,end),evenpay','o');
hold on;
plot(pots(:,end),crosses','ro');
fprintf('Mean Winnigs %.3f +- %.3f',mean(pots(:,end)),std(pots(:,end)));
ttest(pots(:,end))
xlabel('net gain');
ylabel('gain for constant bet');
all.edge = edge;
all.payout = payout;
all.betvals = bets;
all.bets = allbets;
all.meanbets = meanbets;
all.evenpay = evenpay;
all.crosses = crosses;
all.nwins = nwins;

subplot(2,2,2);
hold off;
plot(all.crosses,pots(:,end),'o');
hold on;
plot(all.crosses,all.evenpay,'ro');
xlabel('crossings');
ylabel('winnings, wins');
for j = 1:max(all.crosses)-1;
    x(j) = j+0.5;
    id = find(all.crosses ==j | all.crosses == j+1);
    y(j) = mean(pots(id,end));
    z(j) = mean(all.evenpay(id));
    mb(j) = mean(all.meanbets(id));
    nw(j) = mean(sum(all.nwins(id)));
end
plot(x,y,'b-');
plot(x,z,'r-');
plot(x,mb,'g-');
plot(x,nw,'c-');

all.nbets = [sum((all.bets < bets(2)),2) sum((all.bets > bets(1)),2)];
all.evens = (sum(all.nwins,2) * (all.payout-1) - (runl-sum(all.nwins,2))) * mean(all.meanbets);
subplot(2,2,3);
hold off;
cutoff = std(pots(:,end))./8;
id = find(abs(pots(:,end)) < cutoff);
plot(all.crosses(id)',pots(id,end)-all.evens(id),'o');
%plot(all.crosses,pots(:,end) - all.evens','o');
hold on;
plot(all.crosses(id),pots(id,end) - all.evenpay(id)','r+');
all.bigbets = sum(all.bets' > bets(1));
%plot(all.bigbets',pots(:,end)-all.evens);
xlabel('crossings');
ylabel('Prop-fixed');


[y,x] = smhist(pots(:,end),'smooth',0.3);
[my,mx] = smhist(all.evens,'smooth',0.3);
[ay,ax] = smhist(min(pots),'smooth',0.5);
[by,bx] = smhist(max(pots),'smooth',0.5);
subplot(2,2,4);
hold off;
plot(x,y);
hold on;
plot(mx,my,'r');
plot(ax,ay,'g');
plot(bx,by,'g');
xlabel('gain');
ylabel('frequency');

%
% here is what makes it complicated. Every upward crossing is made by winnning a
% small bet. Every downward crossing is made losing a large bet, so the
% effect of crossings is to distribute the wins unevenly, with a higher
% proportion of wins on the low bets...
id = find(all.nbets(:,1) > 0 & all.nbets(:,2) > 1)
mean(all.nwins(id,:)./all.nbets(id,:))


