function ymatch = MatchInd(x,y, varargin)
%ymatch = MatchInd(x,y)
%for each element of X, find element in y that precedes it
% or 0 if no value of y precedes x(i)
%ymatch = MatchInd(x,y) find element in y that exeeeds x(i)
mode = 'before';
j = 1;
while j <= length(varargin)
    if strncmp(varargin{j},'nearest',5)
        mode = 'nearest';
    elseif strncmp(varargin{j},'after',5)
        mode = 'after';
    end
    j = j+1;
end

dx = bsxfun(@minus,x(:)',y(:));

%dx =
%  x1-y1  x2 -y1  x3-y1 ...
%  x1-y2  x2 -y2  x2-y2 ...
%  x1-y3  x2 -y3  x3-y3 ...

if strcmp(mode,'nearest')
    [a,b] = min(abs(dx));
    ymatch = y(b);
else
    dx = dx >= 0;
    ddx = diff(dx);
    if sum(sum(ddx)) < 0
        if strcmp(mode,'after')
            ddx = [dx(1,:); ddx;];
        else
            ddx = [ddx; -dx(end,:)];
        end
        missed = find(sum(ddx <0) ==0);
        [a,b] = find(ddx<0);
    else
        if strcmp(mode,'after')
            ddx = [ddx; dx(end,:)];
        else
            ddx = [dx(1,:); ddx];
        end
        missed = find(sum(ddx >0) ==0);
        [a,b] = find(ddx);
    end
    ymatch(b) = y(a);
    ymatch(missed) = 0;
end