function PlotFullV(V, varargin)
e = [];
for j = 1:length(varargin)
    if isfield(varargin{j},'Trials')
        E = varargin{j};
    end
    j = j+1;
end
for j = 1:length(V.blkstart)
   
    plot([V.blkstart(j) V.blkstart(j)+V.blklen(j).*V.samper],[1 1],'-');
    hold on;
end

if ~isempty(E)
    if isfield(E.Trials,'Result')
        result = [E.Trials.Result];
    else
        result = ones(1,length(E.Trials,1));
    end
    for j = 1:length(E.Trials)
        h = plot([ E.Trials(j).Start(1) E.Trials(j).Start(1) E.Trials(j).End(1) E.Trials(j).End(1)]./10000,[0 0.9 0.9 0],'r');
        if result(j) ~= 1
            set(h,'color','k');
        end
    end
    set(gca,'ylim',[0 1.1]);
end