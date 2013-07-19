function FullV = FixCoilNoise(FullV, varargin)
period = 4.9448;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'sine',4)
        FindAliasFrequency(FullV);
        return;
    end
    j = j+1;
end
pts = 1:length(FullV.V);

for b = 1:length(FullV.blklen)
    if b > 1
        ts = cumsum(FullV.blklen(1:b))+1;
        ts = ts(end-1);
    else
        ts = 1;
    end
    npts = FullV.blklen(b);
    nb = round(npts/30000);
    bs = round(npts/nb);
    for t = ts:bs:ts+npts-bs
        pts = t:t+bs-1;
        x = 1+round(mod(pts./period,1).*period);
        bins = unique(x);
        for j = 1:length(bins)
            id = find(x == bins(j));
            avg(j) = mean(FullV.V(pts(id)));
        end
        newV(pts) = double(FullV.V(pts))-avg(x);
    end

end

FullV.V = newV;

function FindAliasFrequency(FullV)

fs = 83900:1:84000;
V = double(FullV.V(1:FullV.blklen(1)));
x = [1:length(V)]./30000;
freqs = [1:length(x)].*1/x(end);
for j = 1:length(fs)
    sx = sin(x .* 2 * pi * fs(j));
    ff(j,:) = abs(fft(sx));
    sxs(j,:) = sx(1:100);
end
[a,b] = max(ff,[],2);
plot(fs,freqs(b));
a = abs(fft(V));
%imagesc(ff);
