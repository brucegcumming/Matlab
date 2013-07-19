function resp = EyalCP(varargin)
%resp is mean resp for correct reject, false alarm, hits, misses
ntrials = 1000;
signals = 0.5;
noises = 0;
poolsize = 1;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'noise',4)
        j = j+1;
        noises = varargin{j};
    elseif strncmpi(varargin{j},'poolsize',4)
        j = j+1;
        poolsize = varargin{j};
    elseif strncmpi(varargin{j},'signal',4)
        j = j+1;
        signals = varargin{j};
    elseif strncmpi(varargin{j},'ntrials',4)
        j = j+1;
        ntrials = varargin{j};
    end
j = j+1;
end

if poolsize > 1
    PoolFA(ntrials, noises, signals, poolsize);
    return;
end
for k = 1:length(noises)
    noise = noises(k);
for j = 1:length(signals)
    signal = signals(j);
R(:,1) = randn(ntrials,1);
R(:,2) = randn(ntrials,1)+signal;
N = randn(size(R)).*noise;
D = R+N;
normfactor = std(R(:));
R = R./mean(R(:,2));
crit = signal/(2);
x = -3.:0.1:3.;
dx=0.0;
hit = find(D(:,2) > crit);
miss = find(D(:,2) <= crit);
fa = find(D(:,1) > crit);
cr = find(D(:,1) <= crit);
perf = (length(hit)+length(cr))./length(D(:));
a = hist(R(cr,1),x);
resp(1,j,k) = mean(R(cr,1));
bar(x,a);
hold on;
b = hist(R(fa,1),x,0.5);
resp(2,j,k) = mean(R(fa,1));
bar(x+dx,b,0.5,'r');
c = hist(R(hit,2),x);
resp(3,j,k) = mean(R(hit,2));

bar(x,-c,0.5);
d = hist(R(miss,2),x);
resp(4,j,k) = mean(R(miss,2));
bar(x+dx,-d,0.5,'r');
pc = 100*(length(hit)+length(cr))./length(R(:));
ratios(j,k) = (resp(3,j,k)-resp(2,j,k)+resp(4,j,k)-resp(1,j,k))./(2*signal);
fprintf('Hit-Fa %.2f, Miss-Reject %.2f %.1f correct Ratio %.2f\n',resp(3)-resp(2),resp(4)-resp(1),pc,ratios(j,k))
fprintf('Percent Correct %.2f\n',perf);
resp(5,j,k) = mean(R(:,1));
resp(6,j,k) = mean(R(:,2));
resp(7,j,k) = ratios(j,k);
resp(8,j,k) = (resp(3,j,k)-resp(2,j,k))./signal; %hits - FA. 
resp(9,j,k) = (resp(3,j,k)-resp(1,j,k))./signal; %hits - Correct Rejects. 
resp(10,j,k) = (resp(3,j,k)-resp(4,j,k))./signal; %hits - misses. 
resp(11,j,k) = (resp(2,j,k)-resp(4,j,k)); %FA - misses. 
end
end


function PoolFA(ntrials, noises, signals, poolsize);


k = 1;
for j = 1:length(signals)
    signal = signals(j);
    noise = noises(j);
    ZR = randn(ntrials,poolsize);
    SR = randn(ntrials,poolsize)+signal;
    DZ = mean(ZR,2)+randn(ntrials,1).*noise;
    DS = mean(SR,2)+randn(ntrials,1).*noise;

    ZR = ZR./mean(SR(:));
    SR = SR./mean(SR(:));
    crit = signal/(2);
    x = minmax([minmax(mean(ZR,2)) minmax(mean(SR,2))]);
    x = x(1):(x(2)-x(1))/100:x(2);
    dx=0.0;
    hit = find(DS > crit);
    miss = find(DS <= crit);
    fa = find(DZ > crit);
    cr = find(DZ <= crit);

    y = hist(mean(ZR(cr,:),2),x);
    plot(x,y);
    hold on;
    y = hist(mean(ZR(fa,:),2),x);
    plot(x,y,'r');
    y = hist(mean(SR(hit,:),2),x);
    plot(x,-y);
    y = hist(mean(SR(miss,:),2),x);
    plot(x,-y,'r');
    
    perf = (length(hit)+length(cr))./(2*ntrials);
    resp(1,j,k) = mean(mean(ZR(cr,:)));
    resp(2,j,k) = mean(mean(ZR(fa,:)));
    resp(3,j,k) = mean(mean(SR(hit,:)));
    resp(4,j,k) = mean(mean(SR(miss,:)));

    pc = 100*(length(hit)+length(cr))./(2*ntrials);
    ratios(j,k) = (resp(3,j,k)-resp(2,j,k)+resp(4,j,k)-resp(1,j,k))./(2*signal);
    resp(5,j,k) = mean(ZR(:));
    resp(6,j,k) = mean(SR(:));
    resp(7,j,k) = ratios(j,k);
    resp(8,j,k) = (resp(3,j,k)-resp(2,j,k))./signal; %hits - FA.
    resp(9,j,k) = (resp(3,j,k)-resp(1,j,k))./signal; %hits - Correct Rejects.
    resp(10,j,k) = (resp(3,j,k)-resp(4,j,k))./signal; %hits - misses.
    resp(11,j,k) = (resp(2,j,k)-resp(4,j,k)); %FA - misses.
    fprintf('Hit-Fa %.2f, Miss-Reject %.2f %.1f correct Ratio %.2f\n',resp(3)-resp(2),resp(4)-resp(1),pc,ratios(j,k))
    fprintf('Percent Correct %.2f FA-Miss %.2f\n',perf,resp(11,j,k));
    for p = 1:poolsize
        cps(p,1) = CalcCP(ZR(fa,p),ZR(cr,p));
        cps(p,2) = CalcCP(SR(hit,p),SR(miss,p));
        famiss(p) = mean(ZR(fa,p))-mean(SR(miss,p)); 
    end
    cpa = mean(cps(:));
    fprintf('Mean CP %.3f, FA-Miss %.3f\n',cpa,mean(famiss));
end
