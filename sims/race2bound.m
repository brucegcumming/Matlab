function result = race2bound(varargin)

ntrials = 10;
plottrials = 0;
sigval = 0;
sd = 1;
seed = 0;
crit = [10 10];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'ntrials',5)
        j = j+1;
        ntrials = varargin{j};
    elseif strncmpi(varargin{j},'crit',3)
        j = j+1;
        crit = varargin{j};
    elseif strncmpi(varargin{j},'sigs',3)
        j = j+1;
        sigval = varargin{j};
    elseif strncmpi(varargin{j},'seed',3)
        j = j+1;
        seed = varargin{j};
    elseif strncmpi(varargin{j},'sd',2)
        j = j+1;
        sd = varargin{j};
    end
    j = j+1;
end

if seed
randn('state',seed)
end
duration = 5000;

for nc = 1:size(crit,1)
for ns = 1:length(sigval)

ts = randn(duration,2,ntrials) .*sd;
sig = ones(size(ts)) .* sigval(ns);
sig(:,2,:) = sig(:,2,:) .* -1;
pools = cumsum(ts+sig);
choice = zeros(1,ntrials);
rt = zeros(1,ntrials);
for j = 1:ntrials
aid = find(pools(:,1, j) > crit(nc,1));
bid = find(pools(:,2,j) > crit(nc,2));
if isempty(aid) & isempty(bid)
    choice(j) = -1;
    rt(j) = duration;
elseif isempty(aid) 
    choice(j) = 2;
    rt(j) = bid(1);
elseif isempty(bid) || aid(1) <= bid(1)
    choice(j) = 1;
    rt(j) = aid(1);
elseif bid(1) <= aid(1)
    choice(j) = 2;
    rt(j) = bid(1);
elseif aid(1) <= bid(1)
    choice(j) = 1;
    rt(j) = aid(1);
else
end
sigmax(j) = max(max(pools(:,:,j)));
if plottrials
    plot(squeeze(pools(:,:,j)));
end
end
good = find(ismember(choice,[1 2])); 
  presp(ns,nc) = sum(choice ==2)./sum(ismember(choice,[1 2]));
  rts(ns,nc) = mean(rt);
  rtall{ns,nc} = rt;
  noresp(ns,nc) = sum(choice ==-1);
end
end
   plot(presp,rts,'o-');
   hold on;

%hist(sigmax);
result.presp = presp;
result.sig = sigval;
result.rts = rts;
result.allrts = rtall;
result.noresp = noresp;
result.noise = sd;