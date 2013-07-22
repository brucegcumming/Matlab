function colors = mycolors(varargin)
%
% returns 126 color triplets, designed to be useful for plotting. The
% default color set is meant to work against a white background.
%
% mycolors('mat') returns the color values that matlab used be default when
% plot(... is called with a NxM matrix
% mycolors('spkcolors') returns colors used for spike drawing grey, red, green, blue...


colors = {[1 0 0], [0 0 1], [0 1 0], [1 0 1], [0.9 0.9 0.], [0 1 1],[0.75 0 0.25],...
        [0 0.5 0], [0.5 0 0], [0 0 0.5], [0.5 0.25 0 ], [0.25 0  0.5],[1 0.5 1],...
        [1 1 0.5], [0. 0.5 0], [0.5 0 0.25], [0 0 0.5],[0 0 0],[1 0.5 0.5],[0.5 1. 0.5 ],...
        [0.5 0.5 1],[1 0 0], [0 1 0], [0 0 1], [1 1 0], [1 0 1], [0 1 1],[0.75 0 0.25],[0 0.5 0], ...
        [0.5 0 0], [0 0 0.5], [0.5 0.25 0 ], [0.25 0  0.5],[1 0.5 1],[1 1 0.5], [0.5 0 0], [0 0.5 0],...
        [0 0 0.5],[0 0 0],[1 0.5 0.5],[0.5 1. 0.5 ],[0.5 0.5 1]...
        [1 0 0], [0 1 0], [0 0 1], [1 1 0], [1 0 1], [0 1 1],[0.75 0 0.25],...
        [0 0.5 0], [0.5 0 0], [0 0 0.5], [0.5 0.25 0 ], [0.25 0  0.5],[1 0.5 1],...
        [1 1 0.5], [0.5 0 0], [0 0.5 0], [0 0 0.5],[0 0 0],[1 0.5 0.5],[0.5 1. 0.5 ],...
        [0.5 0.5 1],[1 0 0], [0 1 0], [0 0 1], [1 1 0], [1 0 1], [0 1 1],[0.75 0 0.25],[0 0.5 0], ...
        [0.5 0 0], [0 0 0.5], [0.5 0.25 0 ], [0.25 0  0.5],[1 0.5 1],[1 1 0.5], [0.5 0 0], [0 0.5 0],...
        [0 0 0.5],[0 0 0],[1 0.5 0.5],[0.5 1. 0.5 ],[0.5 0.5 1]...
        [1 0 0], [0 1 0], [0 0 1], [1 1 0], [1 0 1], [0 1 1],[0.75 0 0.25],...
        [0 0.5 0], [0.5 0 0], [0 0 0.5], [0.5 0.25 0 ], [0.25 0  0.5],[1 0.5 1],...
        [1 1 0.5], [0.5 0 0], [0 0.5 0], [0 0 0.5],[0 0 0],[1 0.5 0.5],[0.5 1. 0.5 ],...
        [0.5 0.5 1],[1 0 0], [0 1 0], [0 0 1], [1 1 0], [1 0 1], [0 1 1],[0.75 0 0.25],[0 0.5 0], ...
        [0.5 0 0], [0 0 0.5], [0.5 0.25 0 ], [0.25 0  0.5],[1 0.5 1],[1 1 0.5], [0.5 0 0], [0 0.5 0],...
        [0 0 0.5],[0 0 0],[1 0.5 0.5],[0.5 1. 0.5 ],[0.5 0.5 1]...
    };


spkcolors{1} = [0.5 0.5 0.5];
spkcolors {2} = [1 0 0];
spkcolors {3} = [0 1 0];
spkcolors {4} = [0 0 1];
spkcolors {5} = [1 0 1];
spkcolors {6} = [1 1 0];
spkcolors {7} = [0 1 1];
spkcolors {8} = [0 1 0];
spkcolors {9} = [1 0 0.5];
spkcolors {10} = [0.5 1 0];
spkcolors {11} = [0 0.5 1];
spkcolors {12} = [1 0.5 0.5];
spkcolors {13} = [0.5 1 0.5];
spkcolors {14} = [0.5 0.5 1];
spkcolors {15} = [1 0 1];
spkcolors {16} = [1 1 0];
spkcolors {17} = [0 1 1];
spkcolors {18} = [0 1 0];
spkcolors {19} = [1 0 0.5];
spkcolors {20} = [1 0 0.5];
spkcolors {21} = [1 0 0.5];
spkcolors {22} = [1 0.5 0.5];
spkcolors {23} = [0.5 1 0.5];
spkcolors {24} = [0.5 0.5 1];
spkcolors {25} = [1 0 1];
spkcolors {26} = [1 1 0];
spkcolors {27} = [0 1 1];
spkcolors {28} = [0 1 0];
spkcolors {29} = [1 0 0.5];


colors24 = {[1 0 0], [0 0 1], [0 1 0], [1 0 1], [0.8 0.8 0.0], [0 1 1],[0.75 0 0.25],...
        [0 0.5 0],  [0 0 0.5], [0.5 0.25 0], [0 0 0.5], [0.5 0.5 1],[1 0.5 1],...
        [0.8 0.8 0.5], [0 0.5 0], [0.2 0.4 0.2], [0.2 0.5 1],[0 0 0],[1 0.5 0.5],[0.5 0.5 1 ],...
        [0.5 0.5 1],[1 0 0], [0 1 0], [0 0 1]};

j = 1;
while j <= nargin
    if strncmpi(varargin{j},'mat',3)
        mcolors = [
            0         0    1.0000;
            0    0.5000         0;
            1.0000         0         0;
            0    0.7500    0.7500;
            0.7500         0    0.7500;
            0.7500    0.7500         0;
            0.2500    0.2500    0.2500;
            0         0    1.0000;
            0    0.5000         0;
            1.0000         0         0;
            0    0.7500    0.7500;
            0.7500         0    0.7500;
            0.7500    0.7500         0;
            0.2500    0.2500    0.2500;
            0         0    1.0000;
            0    0.5000         0;
            1.0000         0         0;
            0    0.7500    0.7500;
            0.7500         0    0.7500;
            0.7500    0.7500         0;
            ];
     for k = 1:size(mcolors,2)
         colors{k} = mcolors(k,:);
     end
    elseif strncmpi(varargin{j},'spkcolors',7)
        colors = spkcolors;
    elseif strncmpi(varargin{j},'ncolors',3)
        nc = varargin{j};
        if nc == 24
            colors = colors24;
        end
    elseif strncmpi(varargin{j},'plot',3)
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            nc = varargin{j};
        else
            nc = length(colors)
        end
        if nc == 24
            colors = colors24;
        end
        rows = round(sqrt(nc));
        hold off;
        for k = 1:nc
            plot(rem(k-1,rows),floor((k-1)/rows),'o','color',colors{k},'markerfacecolor',colors{k})
            hold on;
        end
    end
    j = j+1;
end

