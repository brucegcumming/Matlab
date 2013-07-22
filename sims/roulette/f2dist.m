function d = f2dist(f, n, varargin)


xvals = [1:length(f)] - length(f)/2;
j = 1;
while j <= nargin -2
    if strncmpi(varargin{j},'xvals',2)
        j = j+1;
        xvals = varargin{j};
    end
    j = j+1;
end
f = f./sum(f);
nsmp = round((f .* n));
ix = 1;
d = [];
for j = 1:length(f) 
    d(ix:ix+nsmp(j)) = xvals(j); 
    ix = ix+nsmp(j);
end