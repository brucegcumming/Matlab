function [edges, sd, esd, lengths, lengthsd] = RunAngleDists(varargin);

nrep = 10;
sd = pi/2:pi/8:pi;
fita.params = [ 629.2602  257.6254    1.3307  610.0260   53.4942    1.4620];

j = 1;
while j <= nargin
    if strncmpi(varargin{j},'full',3)
        nrep = 10000;
        sd = pi/2:pi/20:3*pi/2;
    elseif strncmpi(varargin{j},'roba',4)
        fita.fitted = FitDist([-pi*4:0.01:pi*4],fita.params,'eval');
    end
    j = j+1;
end

for j = 1:length(sd)
    [a,b,c,d,e] = AngleDist(sd(j),'noplot','nrep',nrep);
    edges(j) = mean(c);
    lengths(j) = mean(e);
    esd(j) = std(c);
    lengthsd(j) = std(e);
end
plot(sd,edges);

    