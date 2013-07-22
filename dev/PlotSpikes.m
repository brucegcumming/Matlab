function PlotSpikes(varargin)


nprobe = 0;
step = 0.5;
j = 1;
while j <= length(varargin)
    if isstruct(varargin{j})
        nprobe = nprobe+1;
        spks{nprobe} = varargin{j};
    elseif strcmpi(varargin{j},'id')
        j = j+1;
        ids{nprobe} = varargin{j}; 
    end
    j = j+1;
end

hold off;
colors = mycolors;
hold off;
for j = 1:length(ids{1})
    for k = 1:nprobe
        plot(spks{k}.values(ids{k}(j),:)+(k-1).*step,'color',colors{k});
        hold on;
    end
end
