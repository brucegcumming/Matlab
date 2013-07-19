function [win, se, edge, angles, radii, ps] = AngleDist(sd, varargin)
% Returns winnings, edge in % for a wrapped Gaussian distribution of
%prediction errors
j = 1;
nloops = 10;
nrep = 1000;
plotdf = 0;
aoffset = 0;
colors = mycolors;
betw = 5;
disttype = 0;
while j < nargin
    if strncmpi(varargin{j},'nloop',3)
        j = j+1;
        nloops = varargin{j};
    elseif strncmpi(varargin{j},'betw',3)
        j = j+1;
        betw = varargin{j};
    elseif strncmpi(varargin{j},'dist',3)
        j = j+1;
        disttype = 1;
        fdist = varargin{j};
        if size(fdist,2) > 1
            xvals = fdist(2,:);
            fdist = fdist(1,:);
        end
    elseif strncmpi(varargin{j},'nrep',3)
        j = j+1;
        nrep = varargin{j};
    elseif strncmpi(varargin{j},'noplot',3)
        plotdf = -1;
    elseif strncmpi(varargin{j},'noplot',3)
        plotdf = -1;
    elseif strncmpi(varargin{j},'type',3)
        j = j+1;
        plotdf = varargin{j};
    elseif strncmpi(varargin{j},'angle',3)
        j = j+1;
        aoffset = varargin{j};
    end
    j = j+1;
end

ps = [];
if plotdf >= 0 
    subplot(2,1,2);
    hold off;
end


if disttype ==1 %%build lookup tablle for uniform rnd -> dist rnd
    nbins = 1000;
    cs = cumsum(fdist);
    [B,I] = unique(cs);
    csi = interp1(cs(I),xvals(I),1:max(cs)/nbins:max(cs));
end

for k = 1:length(aoffset)
for j = 1:length(sd)
    for loop = 1:nloops
        if disttype == 1
            rnd = ceil(rand(nrep,1) .* length(csi));
            a = csi(rnd);
        else
            rnd = randn(nrep,1);
            a = rnd .* sd(j);
        end
        a = atan2(sin(a),cos(a)) + aoffset(k);
        allangles(:,loop) = atan2(sin(a),cos(a));
        angles(j,loop) = atan2(mean(sin(a)),mean(cos(a)));
        radii(j,loop) = abs(i * mean(sin(a))+mean(cos(a)));
        %     p * 36 - (1-p)  = p * 37 -1
        [edge(j,loop), tmp, ps(j,loop)] = CalcEdge(a,'betw',betw);
        if plotdf == 3
            [r, th] = smhist([a a+2*pi],'sd',0.1);
            id = find(th > 0 & th <= 2 * pi);
            [x, y] = pol2cart(th(id),r(id));
            plot(x,y);
            axis image;
        end
    end
    win(j,k) = mean(edge(j,:));
    se(j,k) = std(edge(j,:));
    
    if plotdf >= 0
    subplot(2,1,2);
    if k == 1
        hold off;
    end
    if plotdf == 0
        [r, th] = smhist([angles(j,:) angles(j,:)+2*pi],'sd',0.1,'range',[0:0.01:2*pi]);
    elseif plotdf == 1
        [r, th] = smhist([allangles(:) allangles(:)+2*pi],'sd',0.1,'range',[0:0.1:6.3]);
    end
    %id = find(th >= 0 & th <= 2 * pi);
    [x, y] = pol2cart(th,r);
    plot(x,y,'color',colors{j});
    hold on;
    axis image;
    end
end
rs(:,k) = mean(radii');
rse(:,k) = std(radii');
if plotdf == 3
    subplot(1,1,1);
    for j = 1:length(win);
        plot(radii(j,:),edge(j,:),'o','color',colors{j});
        hold on;
    end
elseif plotdf >= 0
    subplot(2,1,1);
    if k == 1
        hold off;
    end
    errorbar(rs(:,k),win(:,k),se(:,k),'color',colors{k});
    hold on;
end
labels{k} = sprintf('A %.3f',aoffset(k));
end
if plotdf >= 0
    legend(labels);
end