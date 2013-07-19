function phasesf(varargin)

j = 1;
while(j < nargin)
    if(strncmpi(varargin{j},'test',4))
        test = varargin{j};
    end
    j = j+1;
end

phase = [5, 180];
a = [];
b = [];
freqs = [0.2:0.01:10];
for f = freqs
    a = [a phase(1)/( 2 * pi * f)];
    b = [b phase(2)/( 2 * pi * f)];
end

hold off;
plot(freqs,a);
hold on;
plot(freqs,b);
set(gca,'Yscale','Log');

set(gca,'Xscale','Log');
