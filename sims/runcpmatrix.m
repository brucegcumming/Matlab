function [cps, details] = runcpmatrix(varargin)
%run set of simulatons for review figure 2.
%runcpmatrix('runid', x)  runs only congfigurations in the vector x

poolsize = 1000;
n =1;
np = 0;
ncp = poolsize * 2;
runid = [];

j = 1;
while j <= length(varargin)
    if isfield(varargin{j},'cdiff')
        PlotCPRes(varargin{j});
        np = np+1;
    elseif strncmpi(varargin{j},'ncp',3)
        j = j+1;
        ncp = varargin{j};
    elseif strncmpi(varargin{j},'poolsize',5)
        j = j+1;
        poolsize = varargin{j};
    elseif strncmpi(varargin{j},'weights',5)
        [cps, details] = runweights(poolsize, ncp, varargin{:});
        return;
    elseif strncmpi(varargin{j},'runid',5)
        j = j+1;
        runid = varargin{j};
    end
    j =j+1;
end

if np %just plotting
    return;
end


% 25, 27 have systematic difference for non-caual, with
% 27 having lower CP.  26-29 all have CP below the expected line
%
%  a  b  e  f
% c  d  g  h
% e   f  a  b
% g  h  c  d
%              a       b     c       d     e       f       g     h
M(n,:) = [0.05 0.00 0.00 0.05 0.00 0.00 0.00 0.00]; types(n) = 1; n = n+1;
M(n,:) = [0.10 0.00 0.00 0.10 0.00 0.00 0.00 0.00]; types(n) = 1; n = n+1;
M(n,:) = [0.15 0.00 0.00 0.15 0.00 0.00 0.00 0.00]; types(n) = 1; n = n+1;
M(n,:) = [0.05 0.05 0.05 0.05 0.00 0.00 0.00 0.00]; types(n) = 2; n = n+1; %a+b+c
M(n,:) = [0.10 0.10 0.10 0.10 0.00 0.00 0.00 0.00]; types(n) = 2; n = n+1;
M(n,:) = [0.15 0.15 0.15 0.15 0.00 0.00 0.00 0.00]; types(n) = 2; n = n+1;
M(n,:) = [0.15 0.00 0.00 0.15 0.15 0.00 0.00 0.15]; types(n) = 3; n = n+1;%a+d vs e+h
M(n,:) = [0.15 0.00 0.00 0.15 0.05 0.00 0.00 0.05]; types(n) = 3; n = n+1;%a+d vs e+h
M(n,:) = [0.15 0.00 0.00 0.15 0.00 0.15 0.15 0.00]; types(n) = 1; n = n+1;
M(n,:) = [0.10  0.10 0.10 0.00 0.00 0.00 0.00 0.00]; types(n) = 4; n = n+1; % a+b only
M(n,:) = [0.20  0.20 0.20 0.00 0.00 0.00 0.00 0.00]; types(n) = 4; n = n+1; % a+b only
M(n,:) = [0.00  0.10 0.00 0.00 0.10 0.00 0.00 0.00]; types(n) = 5;n = n+1; % a+e
M(n,:) = [0.00  0.10 0.00 0.00 0.00 0.10 0.00 0.00]; types(n) = 5; n = n+1; % a +f 
M(n,:) = [0.00  0.10 0.10 0.00 0.10 0.00 0.00 0.10];types(n) = 6; n = n+1; %b+c and e+h
M(n,:) = [0.00 0.10 0.10 0.00 0.00 0.00 0.00 0.00]; types(n) = 6; n = n+1;  %b+c only
M(n,:) = [0.00 0.20 0.20 0.00 0.00 0.00 0.00 0.00];types(n) = 7; n = n+1;
M(n,:) = [0.10 0.00 0.00 0.05 0.00 0.00 0.00 0.00];types(n) = 7; n = n+1; %a only
M(n,:) = [0.10 0.15 0.15 0.10 0.00 0.00 0.00 0.00]; types(n) = 7; n = n+1; %
M(n,:) = [0.10 0.15 0.15 0.10 0.00 0.00 0.00 0.00]; types(n) = 7; n = n+1;
M(n,:) = [0.00 0.00 0.00 0.10 0.10 0.00 0.00 0.00]; types(n) = 7; n = n+1; %c+d only
M(n,:) = [0.00 0.10 0.10 0.10 0.00 0.00 0.10 0.10]; types(n) = 7; n = n+1;%c+d vs g+h
M(n,:) = [0.20 0.20 0.20 0.20 0.00 0.00 0.00 0.00]; types(n) = 7; n = n+1;
M(n,:) = [0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00]; types(n) = 7; n = n+1;
M(n,:) = [0.05 0.05 0.05 0.05 0.00 0.00 0.00 0.00]; types(n) = 102; n = n+1; %a+b+c
M(n,:) = [0.10 0.10 0.10 0.10 0.00 0.00 0.00 0.00]; types(n) = 102; n = n+1;
M(n,:) = [0.15 0.15 0.15 0.15 0.00 0.00 0.00 0.00]; types(n) = 102; n = n+1;
M(n,:) = [0.15 0.10 0.10 0.15 0.00 0.00 0.00 0.00]; types(n) = 103; n = n+1;
M(n,:) = [0.15 0.05 0.05 0.15 0.00 0.00 0.00 0.00]; types(n) = 103; n = n+1;
M(n,:) = [0.15 0.00 0.00 0.15 0.00 0.00 0.00 0.00]; types(n) = 103; n = n+1;
M(n,:) = [0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00]; types(n) = 102; n = n+1; %a+b+c
details.M = M;
details.type = types;

nloops = round(ncp./(2 * poolsize));
aid = [];
bid = [];
for j = 0:nloops-1
aid = [aid 2*j*poolsize+[1:poolsize/2  poolsize+1:poolsize*1.5]];
end
bid = [aid+poolsize/2];

if isempty(runid)
    runid = 1:size(M,1);
end

for j = runid
    CM = [M(j,1:2) M(j,5:6);
               M(j,3:4) M(j,7:8);
               M(j,5:6) M(j,1:2);
               M(j,7:8) M(j,3:4)];
       if ncp > 0    
       if types(j) > 100 %ignore group 2, as in Cohen+Newsome
           W = ones(1,poolsize);
           W(poolsize/2:end) = 0;
           a = cpsim(poolsize, 'corr', CM, 'weights', W,'ncp',ncp)
       else
           a = cpsim(poolsize, 'corr', CM,'ncp',ncp);
       end
       cps(j,1) = mean(a(aid));
       cps(j,2) = mean(a(bid));
       else
           cps(j,1) = 0;
           cps(j,2) = 0;
       end
       if types(j) > 100
       details.cdiff(j,1) = mean(CM(1,1)) - mean(CM(1,3));
       details.cdiff(j,2) = mean(CM(2,1)) - mean(CM(2,3));
       else
       details.cdiff(j,1) = mean(CM(1,1:2)) - mean(CM(1,3:4));
       details.cdiff(j,2) = mean(CM(2,1:2)) - mean(CM(2,3:4));
       end
end
details.cps = cps;       

function [cps, details] = runweights(poolsize, ncp, varargin)


np = 0;
n = 1;
runid = [];

 j = 1;
while j <= length(varargin)
    if isfield(varargin{j},'cdiff')
        PlotCPRes(varargin{j});
        np = np+1;
    elseif strncmpi(varargin{j},'ncp',3)
        j = j+1;
        ncp = varargin{j};
    elseif strncmpi(varargin{j},'poolsize',5)
        j = j+1;
        poolsize = varargin{j};
    elseif strncmpi(varargin{j},'runid',5)
        j = j+1;
        runid = varargin{j};
    end
    j =j+1;
end

%
%  a  b  e  f
% c  d  g  h
% e   f  a  b
% g  h  c  d
%              a       b     c       d     e       f       g     h
Weights(n,:) = [1 1];
M(n,:) = [0.15 0.00 0.00 0.15 0.00 0.00 0.00 0.00]; types(n) = 2; n = n+1;
Weights(n,:) = [1 1];
M(n,:) = [0.15 0.15 0.15 0.15 0.00 0.00 0.00 0.00]; types(n) = 2; n = n+1;
Weights(n,:) = [1 1];
M(n,:) = [0.075 0.075 0.075 0.075 0.00 0.00 0.00 0.00]; types(n) = 2; n = n+1;

if isempty(runid)
    runid = 1:size(M,1);
end

nloops = round(ncp./(2 * poolsize));
aid = [];
bid = [];
for j = 0:nloops-1
aid = [aid 2*j*poolsize+[1:poolsize/2  poolsize+1:poolsize*1.5]];
end
bid = [aid+poolsize/2];

for j = runid
    CM = [M(j,1:2) M(j,5:6);
               M(j,3:4) M(j,7:8);
               M(j,5:6) M(j,1:2);
               M(j,7:8) M(j,3:4)];
       W = ones(1,poolsize) .* Weights(j,1);
       W(poolsize/2:end) = Weights(j,2);
       a = cpsim(poolsize, 'corr', CM, 'weights', W,'ncp',ncp);
       cps(j,1) = mean(a(aid));
       cps(j,2) = mean(a(bid));
       details.cdiff(j,1) = (sum(CM(1,1:2).*Weights(j,:)) - sum(CM(1,3:4).*Weights(j,:)))./sum(Weights(j,:));
       details.cdiff(j,2) = (sum(CM(2,1:2).*Weights(j,:)) - sum(CM(2,3:4).*Weights(j,:)))./sum(Weights(j,:));
end
details.weights = Weights;
details.M = M;
details.cps = cps;       
details.type = types;


function HitPoint(src, data, id, type)

fprintf('Hit %d,%d\n',id, type);

function PlotCPRes(X, varargin)

hold off;
symbols = '+os^+o^s+o^s';
colors = 'kkkkrrrrbbbb';

rnd = rand(size(X.cps,1),1) .* 0.01;
types = unique(X.type);
for j = 1:length(types)
    id = find(X.type == types(j));
    s = j;
    if ismember(types(j),[102 103])
        s = 7;
    end
    for k = 1:length(id)
        plot(X.cdiff(id(k),1),X.cps(id(k),1),symbols(s),'color',colors(s),'buttondownfcn',{@HitPoint, id(k) , 1});
        hold on;
        if types(j) == 102
            plot(X.cdiff(id(k),2),X.cps(id(k),2),'go','buttondownfcn',{@HitPoint, id(k) , 2});
        elseif types(j) == 103  %don't plot these yet.
        else
            plot(X.cdiff(id(k),2),X.cps(id(k),2),symbols(s),'color',colors(s),'buttondownfcn',{@HitPoint, id(k) , 2});
        end
    end
end

xlabel('mean correlation difference');
ylabel('cp');
figify(gcf,gca,'fontsize',12,'fontname','Ariel');
set(gca,'ylim',[0.48 0.7],'xlim',[-0.01 0.21]);

