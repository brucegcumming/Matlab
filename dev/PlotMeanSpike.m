function PlotMeanSpike(C, varargin)
    addstr = [];
    plots = [ 1 1];
    j = 1; 
    while j <= length(varargin)
        if strncmpi(varargin{j},'addtitile',5)
            j = j+1;
            addstr = varargin{j};
        elseif strncmpi(varargin{j},'lineonly',5)
            plots = [0 1];
        elseif strncmpi(varargin{j},'imageonly',7)
            plots = [1 0];
        end
        j = j+1;
    end
        
    if sum(plots) > 1
    subplot(1,2,1);
    end
    if plots(1)
    hold off;
    imagesc(C.MeanSpike.ms);
    if isfield(C,'probe');
        p = C.probe(1);
    else 
        p = 0;
    end
    line([0 5],[p p],'color','r');
    title(sprintf('P%d Ex %d Gm %.2f (%.2f) %s',p,C.exptno,C.mahal(1),C.mahal(2),addstr));
    if sum(plots) > 1
    subplot(1,2,2);
    end
    end
    if plots(2)
    hold off;
    v = std(C.MeanSpike.ms');
    id = find(v > max(v)/2);
    
    for j = id
        plot(C.MeanSpike.ms(j,:),'r');
        hold on;
        if isfield(C.MeanSpike,'dp') && size(C.MeanSpike.dp,1) >= j
            plot(C.MeanSpike.dp(j,:),'g');
        end
        plot(C.MeanSpike.mu(j,:),'color',[0.5 0.5 0.5]);
    end
    end