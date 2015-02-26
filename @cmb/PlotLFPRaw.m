function PlotLFPRaw(state, Trial, crrate);

times = [1:size(Trial.LFP,1)] .* crrate;
times = times - (Trial.Start(1)-Trial.lfptime)/10000;
plot(times,Trial.LFP);


