function [z, rndamps, meanamp] = nsines(n, varargin)
% [z, rndamps] = nsines(n, varargin)
%
% generates the sum of n sinusoids, frequencies 1...n
%
 % 'amps', x       
%               supplies vector x of amplitudes for each
 % 'phases', x       
%               supplies vector x of phases for each
 % 'freqs', x       
%               supplies vector x of frequencies for each


colors = mycolors;
phases = zeros(n,1);
amps = ones(n,1);
freqs = 1:n;
rndruns = 1;
rndphase = 0;
showplot = 1;
normalise = 1;
rndonoff = 0;
ps = [];

j = 1;
while j < nargin
    if strncmpi(varargin{j},'amps',4)
        j = j+1;
        amps = varargin{j};
    elseif strncmpi(varargin{j},'phases',6)
        j = j+1;
        phases = varargin{j};
    elseif strncmpi(varargin{j},'probs',2)
        j = j+1;
        ps = varargin{j};
        for k = 1:length(ps)
            probs(k,:) = binornd(1,ps(k),rndruns,1);
        end
    elseif strncmpi(varargin{j},'rndphase',4)
        j =j + 1;
        rndphase = 1;
        rndruns = varargin{j};
        if rndruns > 20
            showplot = 0;
        end
    elseif strncmpi(varargin{j},'freqs',5)
        j = j+1;
        freqs = varargin{j};
    elseif strncmpi(varargin{j},'onoff',5)
        rndonoff = 1;
        if (j+1) < nargin & isnumeric(varargin{j+1})
            j = j+1;
            rndonoff = varargin{j};
        end
    elseif strncmpi(varargin{j},'indep',5)
        if rndonoff
            rndonoff = 2;
        end
    elseif strncmpi(varargin{j},'phase',5)
        j = j+1;
        phases = ones(n,1) .* varargin{j};
    elseif strncmpi(varargin{j},'scale',5)
        j = j+1;
        scale =  varargin{j};
    end
    j = j+1;
end

amps = amps * scale;
if isempty(ps)
   probs = ones(length(amps),rndruns); 
end

x = -0.5:0.01:0.5;
for k = 1:rndruns
z = [];
if rndphase
    phases = rand(size(phases)) * pi * 2;
end

if rndonoff == 1
    halfn = floor(length(freqs)/2);
    pf = randperm(length(freqs));
    amps = zeros(size(freqs));
    amps(pf(1:halfn)) = 2 * rndonoff * scale;
    normalise = 0;
elseif rndonoff == 2
    pf = rand(size(freqs));
    amps = (pf > 0.5) * 2 * scale;
    normalise = 0;
end


for j = 1:n
y = probs(j,k) .* amps(j) .* sin(x * 2 * freqs(j) * pi + phases(j));
if showplot
    plot(x,y,'color',colors{j});
   hold on;
end
z = [z;y];
end
if showplot
   plot(x,mean(z),'k');
end

if normalise
    ts(k,:) = mean(z)./mean(amps);
else
    ts(k,:) = mean(z);
end
rndamps(k) = max(abs(ts(k,:)));
meanamp(k) = mean(amps);
end
