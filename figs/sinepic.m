function im = sinepic(sz, varargin)

base = 0.01;
j = 1;
if isnumeric(varargin{1})
    freqs = varargin{1};
end
phases(size(freqs)) = 0;
transparent = 0;
slope = 0;



while j < nargin
    if strncmpi(varargin{j},'freqs',5)
        j = j+1;
        freqs = varargin{j};
        phases = zeros(size(freqs));
    elseif strncmpi(varargin{j},'phases',6)
        j = j+1;
        phases = varargin{j};
    elseif strncmpi(varargin{j},'slope',5)
        j = j+1;
        slope = varargin{j};
    elseif strncmpi(varargin{j},'transparent',6)
        transparent = 1;
        if j +1 < nargin & isnumeric(varargin{j+1})
            j = j+1;
            freqs = varargin{j};
        end
    end
    j = j+1;
end
x = 1:sz;
for j = 1:length(freqs)
    y(j,:) = sin(phases(j) + 2 * pi * freqs(j) .* x./sz);
end
y = mean(y,1);
y = y - min(y);
y = base + (1-base) .* y ./ (range(y));  %% scale 0 ->1;
im = repmat(y,sz,1);

if(slope > 0)
    for j = 1:sz
        left = sz - (j-1) * slope;
        if(left > 0)
            im(j,1:left) = 0;
        end
        right = (sz-1) - (sz-j) * slope;
        if(right > 0)
            im(j,sz-right:sz) = 0;
        end
    end
end

if transparent == 2
   im(:,2:2:end) = 0;
elseif transparent
   im(:,1:2:end) = 0;
end

imagesc(im);
colormap('gray');
