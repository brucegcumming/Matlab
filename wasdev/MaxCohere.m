function coh = MaxCohere(C,varargin)

freqs = 40:60;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'freqs',4)
        j = j+1;
        freqs = varargin{j};
    end
    j = j+1;
end

A = squeeze(sum(C(freqs,:,:)));
[x,y] = meshgrid([1:size(A,1)],1:size(A,2));

subplot(2,1,1);
hold off;
for j = 1:size(A,1)

xi = 1:2*j-1;
yi = 2*j-1:-1:1;

zi = interp2(x,y,A,xi,yi);
%%
plot(xi-j,zi);
hold on;
id = find(abs(xi-j) < 2.5);
coh(j) = sum(zi(id));
end
subplot(2,1,2);
plot(coh);
[a,b] = max(coh);
title(sprintf('Max at %d',b));