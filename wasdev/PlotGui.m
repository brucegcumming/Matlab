function varargout = PlotGui(X,Y, Z, varargin)
%PlotGui(X,Y, Z, varargin)
%Skeleton scirpt that plots X vs Y, then when the mouse is clicked over
%a datapoint, reports the associated Z value (which could be a filename,
%for instance

j = 1;
colorscheme = 'none';
include = ones(size(X));
while j <= length(varargin)
    if strncmpi(varargin{j},'include',6)
        j = j+1;
        include = varargin{j};
    elseif strncmpi(varargin{j},'colors',6)
        j = j+1;
        colorscheme = varargin{j};
    end
    j = j+1;
end
    
    
colors{1} = [0 0 1];
colors{2} = [1 0 0];
colors{3} = [0 1 0];
if isnumeric(colorscheme)
    c = colorscheme;
else
    c = ones(size(X));
end


hold off;
for j = 1:size(X,1)
    if include(j)
        plot(X(j),Y(j),'o','buttondownfcn',{@HitPoint, j},'color',colors{c(j)});
        hold on;
    end
end
set(gcf,'UserData',Z);


function HitPoint(a, b, pt)

h = get(a,'parent');
htype = get(h,'type');
while ~strcmp(htype,'figure')
h = get(h,'parent');
htype = get(h,'type');
end
Z = get(h,'UserData');
fprintf('you hit point %d Z = %.1f\n',pt,Z(pt));




