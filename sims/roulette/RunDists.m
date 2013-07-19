function [e,r,p] = RunDists(sd, varargin)
%
%  look at effect of sd/Width with non-random calculation

disttype = 0;
j = 1;
while j <= nargin -1
    if strncmpi(varargin{j},'dist',3)
        j = j+1;
        disttype = 1;
        fdist = varargin{j};
        if size(fdist,2) > 1
            xvals = fdist(2,:);
            fdist = fdist(1,:);
        end
    end
    j = j+1;
end


x = -4*pi:0.001:4*pi;

if disttype == 0
    for j = 1:length(sd)
        y = Gauss(sd(j),x);
        [e(j), r(j), p(j)] = CalcDistEdge(x,y);
    end
elseif disttype == 1;
    for j = 1:length(sd)        
        y = interp1(xvals .*sd(j),fdist,x);
        y(find(isnan(y))) = 0;
        [e(j), r(j), p(j)] = CalcDistEdge(x,y);
    end    
end