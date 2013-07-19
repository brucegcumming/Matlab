function z = NsineSpec(nruns,varargin)
%
% Builds sum of sinewaves with different spectral properties
% to compare different methods for testing adpative filtering

figlabel{1} = 'NsinePlot';
freqs = 1:18;
amps = 1./freqs;
lfreqs = exp(0:log(18)/17:log(18));
scale = 1;
indep = 0;
argon = {};

j = 1;
while j < nargin
    if strncmpi(varargin{j},'indep',3)
        argon = {argon{:} 'indep'};
        indep = 1;
    elseif strncmpi(varargin{j},'freqs',3)
        j = j+1;
        freqs = varargin{j};
        nf = length(freqs);
        amps = 1./freqs;
        lfreqs = exp(0:log(nf)/(nf-1):log(nf));
    elseif strncmpi(varargin{j},'scale',3)
        j = j+1;
        scale = varargin{j};
    end
    j = j+1;
end

GetFigure(figlabel{1});
hold off;
[a,b, c] = nsines(length(freqs),'freqs',freqs,'rndphase',nruns,'noplot','onoff','scale',scale,argon{:});
z(1,:) = b;
[a,b] = hist(b);
bar(b,a,0.5);
[a,b,c] = nsines(length(freqs),'freqs',lfreqs,'rndphase',nruns,'noplot','onoff','scale',scale,argon{:});
hold on;
z(2,:) = b;
[a,b] = hist(b);
bar(b+0.02,a,0.5,'r');
if 0
    [a,b] = nsines(length(freqs),'freqs',freqs,'amps',amps,'rndphase',nruns,'noplot');
    hold on;
    [a,b] = hist(b);
    bar(b+0.04,a,'g');
end
