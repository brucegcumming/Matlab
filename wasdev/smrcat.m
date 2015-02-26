function smrcat(outf, varargin)

nin = 1;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'lfponly',4)
    else
        infiles{nin} = varargin{j};
        nin =nin+1;
    end
    j = j+1;
end

cHeader.files = infiles;
cHeader.TimeOffsets(1) = 0;

LFP = [];
EM = [];
load(infiles{1});
Key = Ch31;
Text = Ch30;
Spikes = Ch5;
StimOn = Ch8;
Frames = Ch7;

lfpfile = strrep(infiles{1},'.mat','.lfp.mat');
if exist(lfpfile)
    load(lfpfile);
    LFPall = LFP;
    if ~isfield(LFP,'Pulses');
        LFPall.Pulses = [];
    end
end

emfile = strrep(infiles{1},'.mat','.em.mat');
if exist(emfile)
    load(emfile);
    EM = Expt;
end

if exist('Ch9','var')
    StimChange = C9;
end

for j = 2:length(infiles)
    toff = max([max(Key.times) max(Text.times) max(Spikes.times)]);
    cHeader.TimeOffsets(j) = toff;

    clear Ch*
    load(infiles{j});
    Text.times = [Text.times; Ch30.times+toff];
    Text.text = [Text.text; Ch30.text];
    Text.codes = [Text.codes; Ch30.codes];
    Key.times = [Key.times; Ch31.times+toff];
    Key.codes = [Key.codes; Ch31.codes];
    Spikes.times = [Spikes.times; Ch5.times + toff];
    Spikes.codes = [Spikes.codes; Ch5.codes];
    Spikes.values = [Spikes.values; Ch5.values];
    StimOn.times = [StimOn.times; Ch8.times + toff];
    StimOn.level = [StimOn.level; Ch8.level];
    Frames.times = [Frames.times; Ch7.times];

    lfpfile = strrep(infiles{j},'.mat','.lfp.mat');
    if exist(lfpfile)
        load(lfpfile);
        for k = 1:length(LFP.Trials)
            LFP.Trials(k).Start = LFP.Trials(k).Start + toff;
            LFP.Trials(k).ftime = LFP.Trials(k).ftime + toff;
        end
        LFPall.Trials = [LFPall.Trials LFP.Trials];
        for k = 1:length(LFP.Trials)
            LFP.Trials(k).Start = LFP.Trials(k).Start + toff;
            LFP.Trials(k).ftime = LFP.Trials(k).ftime + toff;
        end
        if isfield(LFP,'Pulses')
            for k = 1:length(LFP.Trials)
                LFP.Pulses(k).ftime = LFP.Pulses(k).ftime + toff;
            end
        end
    end

    emfile = strrep(infiles{j},'.mat','.em.mat');
    if exist(emfile)
        load(emfile);
        for k = 1:length(Expt.Trials)
            Expt.Trials(k).Start = Expt.Trials(k).Start + toff;
            Expt.Trials(k).ftime = Expt.Trials(k).ftime + toff;
        end
        EM.Trials = [EM.Trials Expt.Trials];
    end
end
    
Ch31 = Key;
Ch30 = Text;
Ch5 = Spikes;
Ch8 = StimOn;
Ch7 = Frames;
save(outf,'Ch5','Ch8','Ch30','Ch31','Ch7','cHeader');

if isempty(strfind(outf,'.mat'))
    outf = [outf '.mat'];
end

if ~isempty(LFPall)
    if isempty(LFPall.Pulses)
        LFP = rmfield(LFPall.Pulses);
    else
        LFP = LFPall;
    end
    save(strrep(outf,'.mat','.lfp.mat'),'LFP','cHeader');
end

if ~isempty(EM)
    Expt = EM;
    save(strrep(outf,'.mat','.em.mat'),'Expt','cHeader');
end


