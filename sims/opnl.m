function opnl(varargin);
%Quick look at how amplitude of generating input (+-A)
%affects amplitude of response with an output nonlinearity
%
%amplidude of response relative to baseline matters
Amps = [1 2 4 8 16 32 64 128];
noises = [0 1 2 4 8];
U = 100;
nr = 100;
G = ones(2,nr);
G(1,:) = -1;
rnd = rand(size(G));
for j = 1:length(Amps)
    a = Amps(j);
for k = 1:length(noises)
    noise = noises(k);
    R = Amps(j) .* G + (rnd.*noise);
    R = (R+U).^2;
    uR = (U + (rnd.*noise)).^2;
    A = mean(R') - mean(uR(:));
    ratio(j,k) = -A(1)/A(2);
    resps(j,k,:) = A;
end
end
plot(Amps./U, ratio);
xlabel('reponse amplitude (ratio of Uncorr)');
ylabel('Response -A/+A');