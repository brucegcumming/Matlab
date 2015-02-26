function [aid, bid] = TwoProbes(sa,sb, varargin)
%
% [aid, bid] = TwoProbes(sa,sb, varargin)
% find simultaneous events on both channels and compare to
% independent picks from each
% ...,'dt',t)  set the min interval for 'simultaneous
% ..,'hists')  plot histograms of height for the sum. 
n = 1;
start = 10000;
ns = 50000;
td = 0.0005;
plotcc = 0;
plotv = 0;
plottype = 0;
dopca = 0;
sign = -1;
PLOTHISTS = 4;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'cc',2)
        plotcc = 1;
    elseif strncmpi(varargin{j},'dt',2)
        j = j+1;
        td = varargin{j};        
    elseif strncmpi(varargin{j},'hists',4)
        plottype = PLOTHISTS;
    elseif strncmpi(varargin{j},'pca',3)
        dopca = 1;
    elseif strncmpi(varargin{j},'sign',4)
        j = j+1;
        sign = varargin{j};
    elseif strncmpi(varargin{j},'start',4)
        j = j+1;
        start = varargin{j};
    elseif strncmpi(varargin{j},'nspkies',3)
        j = j+1;
        ns = varargin{j};
    end
    j = j+1;
end

if sign < 0
    id = find(sa.values(:,9) < 0);
    sa.times = sa.times(id);
    sa.values = sa.values(id,:);
    id = find(sb.values(:,9) < 0);
    sb.times = sb.times(id);
    sb.values = sb.values(id,:);
elseif sign > 0
    id = find(sa.values(:,9) > 0);
    sa.times = sa.times(id);
    sa.values = sa.values(id,:);
    id = find(sb.values(:,9) > 0);
    sb.times = sb.times(id);
    sb.values = sb.values(id,:);
end


[t, id] = min(abs(sa.times(1+start) - sb.times));
if t < td 
    aid(n) = start+1;
    bid(n) = id;
    diffs(n) = t;
    n = n+1;
end

h = zeros(1,201);
for j = 1:ns
    [t, a] = min(abs(sa.times(j+start) - sb.times(id:id+100)));
    id = id+a-1;
    if plotcc
        h = h + hist(sa.times(j+start) -sb.times(id-100:id+100),[-0.05:0.0005:0.05]);
    end
     if t < td
        aid(n) = j+start;
        bid(n) = id;
        diffs(n) = t;
        n = n+1;
    end
end


xc = sum(sa.values(aid,:) .* sb.values(bid,:),2);
%xc = range(sa.values(aid,:)+sb.values(bid,:),2);
subplot(2,1,1);
hist(xc,100);
if plotcc
    bar(h((2:200)));
end
id = find(xc > 8);
%plot(sa.values(aid(id),:)');
subplot(2,1,2);
hold off;
if dopca
    [E,V] = eig(cov([sa.values(aid,:) sb.values(bid,:)]));
    pc(1,:) = (E(1:32,end)'-mean(E(1:32,end))) * sa.values(aid,:)';
    pc(2,:) = (E(33:end,end)'-mean(E(33:end,end))) * sb.values(bid,:)';    
    pc(3,:) = E(1:32,63)' * sa.values(aid,:)';
    pc(4,:) = E(33:end,63)' * sb.values(bid,:)';    
    plot(pc(1,:),pc(2,:),'.');
    hold on;
    plot(pc(1,:)+pc(3,:),pc(2,:)+pc(4,:),'r.');
    a = find(pc(1,:) > 0 & pc(2,:) > 0);
    subplot(2,1,1);
    hold off;
    plot(sa.values(aid(a),:)'-repmat(mean(sa.values(aid(a),:)'),32,1));
    subplot(2,1,2);
    hold off;
    plot(sb.values(bid(a),:)'-repmat(mean(sb.values(bid(a),:)'),32,1));
end
if plotv
plot(sb.values(bid(id),:)');
%plot(sa.values(aid(id),:)',sb.values(bid(id),:)');
else
id = start:start+length(aid);

%heights spikes alone, indep picks
aah = max(sa.values(id,:)') - min(sa.values(id,:)');
bbh = max(sb.values(id,:)') - min(sb.values(id,:)');

%hieght for sum
ph = max((sb.values(bid,:)+sa.values(aid,:))'/2) - ...
    min((sb.values(bid,:)+sa.values(aid,:))'/2);
plot(aah,bbh,'.','markersize',1);
hold on;

%heights for matched spikes alone
ah = max(sa.values(aid,:)') - min(sa.values(aid,:)');
bh = max(sb.values(bid,:)') - min(sb.values(bid,:)');
aic = start:start+length(aid);
bic = start:start+length(bid);

plot(ah,ph,'r.','markersize',2);
if plottype == PLOTHISTS
    [y,x] = smhist(ph);
    hold off;
    plot(x,y);
    hold on;
    [y,x] = smhist(aah);
    plot(x,y,'r');
    [y,x] = smhist(bbh);
    plot(x,y,'g');
    legend('Sum',sa.title,sb.title)
end    
end
title(sprintf('%d/%d spikes',n,ns));