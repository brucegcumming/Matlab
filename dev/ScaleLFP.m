function LFP = ScaleLFP(LFP, varargin)

plottype = 1;
nfreq = 50;
checkdim = 0;
minlen = 0;
chscale = 0;
needft = 0;
fixch = 20;
choff = 0;
checkscale = 0;
exptype = [];


usetrials = [1:length(LFP.Trials)];

j = 1;

while j <= length(varargin)
    if strncmpi(varargin{j},'image',4)
        plottype = 2;
    elseif strncmpi(varargin{j},'checkscale',7)
        checkscale =1;
    elseif strncmpi(varargin{j},'check',4)
        checkdim = 1;
        if length(varargin) > j & isnumeric(varargin{j+1})
            j = j+1;
            minlen = varargin{j};
        end
    elseif strncmpi(varargin{j},'scale',4)
        j = j+1;
        chscale = varargin{j}(1);
        if length(varargin{j}) > 1
        choff = varargin{j}(2);
        end
    elseif strncmpi(varargin{j},'trials',4)
        j = j+1;
        usetrials = varargin{j};
    end
    j = j+1;
end

if checkscale
    hold off;
    for j = usetrials;
        if size(LFP.Trials(j).LFP,2) > 15
        plot(LFP.Trials(j).LFP(:,15),LFP.Trials(j).LFP(:,17),'.');
      hold on;
        end
    end
    return;
end
if ~isfield(LFP.Trials,'co') & strmatch(exptype,'co')
    [Expt.Trials.co] = deal(-1);
end

for j = 1:length(LFP.Trials)
    chs(j) = size(LFP.Trials(j).LFP,2);
    lens(j) = size(LFP.Trials(j).LFP,1);
    if chs(j) > fixch & fixch > 0
        LFP.Trials(j).LFP(:,fixch) = mean(LFP.Trials(j).LFP(:,[fixch-1 fixch+1]),2);
    end
    if chscale > 0 & chs(j) > 16
        LFP.Trials(j).LFP(:,17:chs(j)) = LFP.Trials(j).LFP(:,17:chs(j)) .* chscale + choff; 
    end
end

    