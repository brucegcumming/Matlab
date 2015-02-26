function [x, score] = FindEllipse(xy,idlist, varargin)
%[x, score] = FindEllipse(xy,idlist, varargin) find the ellipse that encloses points
%in xy where idlist ==2, and
%excludes points hwere idlist == 1
guess = [];
plottype = 0;
trackfit = 0;
j = 1;
cluster = 2;
while j <= length(varargin)
    if strncmpi(varargin{j},'guess',4)
        j = j+1;
        guess = varargin{j};
    elseif strncmpi(varargin{j},'cluster',4)
        j = j+1;
        cluster = varargin{j}+1;
    elseif strncmpi(varargin{j},'plot',4)
        plottype =1;
    elseif strncmpi(varargin{j},'tag',3)
        j = j+1;
        GetFigure(varargin{j});
    end
    j = j+1;
end
if diff(size(idlist)) > 0
    idlist = idlist';
end
id = find(idlist == cluster);
nid = find(idlist ~= cluster);

idlist(id)=1;
idlist(nid) = 0;
if isempty(guess)
    xyr([1 2]) = mean(xy(id,:));
    xyr([3 4]) = std(xy(id,:));
    xyr(5) = 0; %angle;
else
    xyr = guess;
end
if trackfit && isappdata(0,'scores')
    rmappdata(0,'scores');
end
[~, score] = TryEllipse(xyr, xy, idlist, trackfit);
options = optimset('MaxFunEvals',100000,'maxiter',1000,'display','off');
if trackfit
    setappdata(0,'scores',[xyr score]);
end
[x, score, exitflag] = fminsearch(@TryEllipse, xyr,options,  xy, idlist, trackfit);
[~, score] = TryEllipse(x, xy, idlist);
scores = [];
diffs = -0.01:0.001:+0.01;
for a = diffs(:)'
    [~, b] = TryEllipse([x(1:4) x(5) + a], xy, idlist, trackfit);
    scores(end+1,:) = b;
end
[a,b] = min(sum(scores,2));
score = scores(b,:);
x(5) = x(5) + diffs(b);
x(5) = x(5).*10;
r = CalcClusterDistance(x,xy);
if plottype
    hold off;
    PlotND(xy,[],'idlist',idlist);
    hold on;
    DrawEllipse(x);
end

function [cost, score] = TryEllipse(xyr, xy, idlist, trackfit)

id = find(idlist ==2);
nid = find(idlist ~= 2);
if nargin < 4
    trackfit = 0;
end

r = CalcClusterDistance(xyr,xy);
score(1) = sum(r <= 1 & ~idlist);
score(2) = sum(r > 1 & idlist);
cost = sum(score);
if trackfit
scores = getappdata(0,'scores');
if ~isempty(scores)
    scores = [scores; xyr score];
    setappdata(0,'scores',scores);
end
end


function r = CalcClusterDistance(xyr, xy)
  xys = xyrotate(xy(:,1)-xyr(1),xy(:,2)-xyr(2),xyr(5).*10);
  r = sqrt(((xys(:,1))./xyr(3)).^2 + ((xys(:,2))./xyr(4)).^2);
