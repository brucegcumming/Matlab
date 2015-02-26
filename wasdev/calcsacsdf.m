function sdf = calcsacsdf(Trials, width, times,flag)
% calcsdf takes a vector list of spike times, spikes
% and calculates a spike density function with either 
% a Gaussian of S.D. width  
% or an exponential of time constant width (if flag = 'exp')
% in steps of step (timestamp units).
%

spikes = [];
nsac = 0;
for j = 1:length([Trials.Trial])
  for sac = [Trials(j).Saccades.start]
    spikes = [spikes Trials(j).Spikes'-sac];
    nsac = nsac + 1;
  end
end

rates = [];
if(strmatch(flag,'exp'));
for t = [times]
  rates = [rates sum(exp((spikes(find(spikes <t)) -t)/width))];
end
sdf = [rates .* 10000/width];
else
for t = [times]
  rates = [rates sum(exp(-((spikes-t).^2)/(2 * width^2)))];
end
sdf = [rates .* 10000/(sqrt(2 * pi) * width * nsac)];
end

    
  
  