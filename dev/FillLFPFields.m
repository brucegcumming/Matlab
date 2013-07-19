function LFP = FillLFPFields(LFP,varargin)
chscale = 20;
if ~isfield(LFP.Stimvals,'e3')
    LFP.Stimvals.e3 = 'e0';
end

if ~isfield(LFP.Stimvals,'Nf')
    LFP.Stimvals.Nf = '120';
end

if ~isfield(LFP.Header,'Name')
    LFP.Header.Name = 'OTRC';
end
if ~isfield(LFP.Trials,'TrialStart')
    for j = 1:length(LFP.Trials)
        LFP.Trials(j).TrialStart = LFP.Trials(j).Start(1);
    end
end

if ~isfield(LFP.Header,'chscale')
    for j = 1:length(LFP.Trials)
        LFP.Trials(j).LFP(:,17:24) = LFP.Trials(j).LFP(:,17:24)./chscale;
    end
    LFP.Header.chscale = chscale;
end
