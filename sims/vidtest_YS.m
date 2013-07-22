function ds = vidtest(varargin);

% simulation of spike syncrhony measures used in Saalmann et al, for
% two independent spike trians (not poision)
% change nspks to simulate different changes in rate for MT/LIP
% vidtest('mtonly')   increases rate only in MT
% vidtest('same')   increases rate in MT and LIP the same
% vidtest('mtmore')   increases rate ib both, but change is larger in MT


triallen = 2.0; %Arbitray scaling of time interval.
%
% nspks(1,1) is # of spikes in MT for contidion 1
% nspkes(1,2) is # of spikes in MT for contidion 2
% nspks(2,1) is # of spikes in LIP for contidion 1
% nspkes(2,2) is # of spikes in LIP for contidion 2


nspks = [50 80; 50 80]; %60% increase in MT and LIP

nloops = 10000;

%delta is counting interval, in units of fraction of triallen
delta = 0.0075;
deltas = [0.0025:0.00125:0.0075];

j = 1;
while j <= nargin
    if strncmpi(varargin{j},'mtonly',3)
        nspks = [50 55; 35 35]; %10% increase in MT only
    elseif strncmpi(varargin{j},'decline',7)
        nspks = [50 41; 35 29]; %18% decrease in MT, 17% in LIP
    elseif strncmpi(varargin{j},'lponly',3)
        nspks = [50 50; 35 39]; %11% increase in LIP only
    elseif strncmpi(varargin{j},'lpmore1',7)
        nspks = [50 55; 35 39]; %10% increase in MT, 11% in LIP
    elseif strncmpi(varargin{j},'lpmore2',7)
        nspks = [50 57; 35 50]; %14% increase in MT, 43% in LIP
    end
    j = j+1;
end


for nd = 1:length(deltas)
    delta = deltas(nd);
    for k = 1:size(nspks,2);
        nspk = nspks(:,k);
        for l = 1:nloops
% distribute spikes at random through the trial. lp and mt are
% vectors of spike times, in units of fraction of triallen.
            lp = rand(nspk(2),1);
            mt = rand(nspk(1),1);
            
            both = 0;
%for each MT spike, count the number of LIP spikes preceding in
% time window delta
            for j = 1:nspk(1)
                pre = mt(j) - lp;
                n = sum(pre > 0 & pre < delta);
                if n %add 1 to coincidence count, MT precede by LIP
                    both = both +1;
                end
            end
 %Now count MT spikes following an LIP spike           
            mtlip(l,k) = both;
            both = 0;
            for j = 1:nspk(2)
                post = sum(mt > lp(j) & mt < lp(j) +delta);
                if post
                    both = both+1; %LIP spike followed by MT
                end
            end
            lipmt(l,k) = both;
        end
    end
    MTpreceded(nd,:) = 100 .* mean(mtlip)./nspks(1,:);
    LIPfollowed(nd,:) = 100 .* mean(lipmt)./nspks(2,:);
end

subplot(2,1,1);
plot(1000 .* triallen .* deltas,MTpreceded,'o-');
xlabel('Interval (ms)');
ylabel('percent(spikes MT) preceded by spike LIP');
nd = nd;
legend({sprintf('%d spks/sec',nspks(1,1)/triallen) sprintf('%d spks/sec',nspks(1,2)/triallen)});
subplot(2,1,2);
plot(1000 .* triallen .* deltas,LIPfollowed,'o-');
ylabel('LIP followed by MT');
legend({sprintf('%d spks/sec',nspks(2,1)/triallen) sprintf('%d spks/sec',nspks(2,2)/triallen)});

