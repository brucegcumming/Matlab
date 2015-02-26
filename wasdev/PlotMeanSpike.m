function PlotMeanSpike(C, varargin)
    addstr = [];
    plots = [ 1 1];
    smoothw = 0;
    j = 1; 
    while j <= length(varargin)
        if strncmpi(varargin{j},'addtitile',5)
            j = j+1;
            addstr = varargin{j};
        elseif strncmpi(varargin{j},'lineonly',5)
            plots = [0 1];
        elseif strncmpi(varargin{j},'imageonly',7)
            plots = [1 0];
        elseif strncmpi(varargin{j},'csdonly',7)
            plots = [2 0];
        elseif strncmpi(varargin{j},'withcsd',7)
            plots = [1 2];
        elseif strncmpi(varargin{j},'withdy',6)
            plots = [1 3];
        elseif strncmpi(varargin{j},'smooth',5)
            smoothw = 1;
        end
        j = j+1;
    end
        
    if sum(plots>0) > 1
        subplot(1,2,1);
    end
    if sum(plots ==2) 
        CSD = diff(diff(C.MeanSpike.ms));
        if smoothw > 0
            [~,~,Z] = gauss2d(smoothw,-5:5);
            CSD = conv2(CSD,Z,'same');
        end
    end
    if sum(plots ==3) 
        dVdY = diff(C.MeanSpike.ms);
        if smoothw > 0
            [~,~,Z] = gauss2d(smoothw,-5:5);
            dVdY = conv2(dVdY,Z,'same');
        end
    end
    if plots(1)
        hold off;
        if plots(1) == 2 %CSD
            Im = CSD;
        else
            Im = C.MeanSpike.ms;
        end
        imagesc(Im);
        if isfield(C,'probe');
            p = C.probe(1);
        else
            p = 0;
        end
    line([0 5],[p p],'color','r');
    title(sprintf('P%d Ex %d Gm %.2f (%.2f) %s',p,C.exptno,C.mahal(1),C.mahal(2),addstr));
    if sum(plots>0) > 1
    subplot(1,2,2);
    end
    end
    if plots(2)
    hold off;
    if plots(2) == 2
        imagesc(CSD);
    elseif plots(2) == 3
        imagesc(dVdY);        
    else
        v = std(C.MeanSpike.ms');
        id = find(v > max(v)/2);
        
        for j = id
            plot(C.MeanSpike.ms(j,:),'r');
            hold on;
            if isfield(C.MeanSpike,'dp') && size(C.MeanSpike.dp,1) >= j
                plot(C.MeanSpike.dp(j,:),'g');
            end
            if isfield(C.MeanSpike,'mu')
                plot(C.MeanSpike.mu(j,:),'color',[0.5 0.5 0.5]);
            end
        end
    end
    end
