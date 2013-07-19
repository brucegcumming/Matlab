function f = loxsim(varargin)
%
% simulate pdfs for color space spanned by the lox/cre recombinase trick

ngenes = 3;
ncolors = 4;
loops = 1;
nsim = 1000000;
loops = 1;
j = 1;
while j <= length(varargin)
     if strncmpi(varargin{j},'loops',2)
        j = j+1;
        loops = varargin{j};
     elseif strncmpi(varargin{j},'ncolors',2)
        j = j+1;
        ncolors = varargin{j};
    elseif strncmpi(varargin{j},'ngenes',2)
        j = j+1;
        ngenes = varargin{j};
    elseif strncmpi(varargin{j},'nsim',2)
        j = j+1;
        nsim = varargin{j};
    end
    j = j+1;
end
allvals = [];
for j = 1:loops
    triplets = floor(rand(nsim,ngenes) .* ncolors);
    vals = sum((ngenes+1).^triplets,2);
    allvals = [allvals; vals];
end
bins = hist(allvals,max(allvals));
id = find(bins > 0);
f = sort(bins(id)./(nsim*loops),'descend');
hold off;
for j = 1:length(f)
h(j) = bar(j,f(j),1);
hold on;
end
if ngenes == 8
set(h(1),'facecolor',[2/8 3/8 3/8]);
set(h(2),'facecolor',[3/8 3/8 2/8]);
set(h(3),'facecolor',[3/8 2/8 3/8]);

set(h(4),'facecolor',[4/8 2/8 2/8]);
set(h(5),'facecolor',[2/8 4/8 2/8]);
set(h(6),'facecolor',[2/8 2/8 4/8]);
j = 7;
set(h(j),'facecolor',[1/8 3/8 4/8]); j = j+1;
set(h(j),'facecolor',[3/8 1/8 4/8]); j = j+1;
set(h(j),'facecolor',[1/8 4/8 3/8]); j = j+1;
set(h(j),'facecolor',[3/8 4/8 1/8]); j = j+1;
set(h(j),'facecolor',[4/8 1/8 3/8]); j = j+1;
set(h(j),'facecolor',[4/8 3/8 1/8]); j = j+1;

set(h(j),'facecolor',[4/8 4/8 0/8]); j = j+1;
set(h(j),'facecolor',[0/8 4/8 4/8]); j = j+1;
set(h(j),'facecolor',[4/8 0/8 4/8]); j = j+1;

set(h(j),'facecolor',[1/8 5/8 2/8]); j = j+1;
set(h(j),'facecolor',[2/8 5/8 1/8]); j = j+1;
set(h(j),'facecolor',[5/8 1/8 2/8]); j = j+1;
set(h(j),'facecolor',[5/8 2/8 1/8]); j = j+1;
set(h(j),'facecolor',[1/8 2/8 5/8]); j = j+1;
set(h(j),'facecolor',[2/8 1/8 5/8]); j = j+1;

set(h(j),'facecolor',[0/8 5/8 3/8]); j = j+1;
set(h(j),'facecolor',[3/8 5/8 0/8]); j = j+1;
set(h(j),'facecolor',[5/8 0/8 3/8]); j = j+1;
set(h(j),'facecolor',[5/8 3/8 0/8]); j = j+1;
set(h(j),'facecolor',[0/8 3/8 5/8]); j = j+1;
set(h(j),'facecolor',[3/8 0/8 5/8]); j = j+1;

set(h(j),'facecolor',[1/8 6/8 1/8]); j = j+1;
set(h(j),'facecolor',[6/8 1/8 1/8]); j = j+1;
set(h(j),'facecolor',[1/8 1/8 6/8]); j = j+1;

set(h(j),'facecolor',[6/8 2/8 0/8]); j = j+1;
set(h(j),'facecolor',[6/8 0/8 2/8]); j = j+1;
set(h(j),'facecolor',[2/8 6/8 0/8]); j = j+1;
set(h(j),'facecolor',[0/8 6/8 2/8]); j = j+1;
set(h(j),'facecolor',[2/8 0/8 6/8]); j = j+1;
set(h(j),'facecolor',[0/8 2/8 6/8]); j = j+1;

set(h(j),'facecolor',[7/8 1/8 0/8]); j = j+1;
set(h(j),'facecolor',[7/8 0/8 1/8]); j = j+1;
set(h(j),'facecolor',[1/8 7/8 0/8]); j = j+1;
set(h(j),'facecolor',[0/8 7/8 1/8]); j = j+1;
set(h(j),'facecolor',[0/8 1/8 7/8]); j = j+1;
set(h(j),'facecolor',[1/8 0/8 7/8]); j = j+1;

set(h(j),'facecolor',[8/8 0/8 0/8]); j = j+1;
set(h(j),'facecolor',[0/8 8/8 0/8]); j = j+1;
set(h(j),'facecolor',[0/8 0/8 8/8]); j = j+1;

end
mean(triplets);
title(sprintf('%d colors (%d)',length(id),max(allvals)));
ylabel('Frequency');
xlabel('Color combination');

