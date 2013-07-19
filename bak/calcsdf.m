function sdf = calcsdf(Trials, width, step,flag)
% calcsdf takes a vector list of spike times, spikes
% and calculates a spike density function with either 
% a Gaussian of S.D. width  
% or an exponential of time constant width (if flag = 'exp')
% in steps of step (timestamp units).
%

spikes = [];
for j = 1:length([Trials.Trial])
  spikes = [spikes Trials(j).Spikes'];
end

rates = [];
times = [];
if(strmatch(flag,'exp'));
for t = -300:step:max(spikes)+width
  rates = [rates sum(exp((spikes(find(spikes <t)) -t)/width))];
  times = [times t];
end
sdf = [times; rates .* 10000/width];
else
for t = -300:step:max(spikes)+width
  rates = [rates sum(exp(-((spikes-t).^2)/(2 * width^2)))];
  times = [times t];
end
sdf = [times; rates .* 10000/(sqrt(2 * pi) * width)];
end

function sdf = circularsdf(spikes, period, width, step,flag)

rates = [];
times = [];
if(strmatch(flag,'exp'));
for t = 0:step:period
  rates = [rates sum(exp((spikes(find(spikes <t)) -t)/width))];
  times = [times t];
end
sdf = [times; rates .* 10000/width];
else
for t = 0:step:period
  rates = [rates sum(exp(-(min([abs(spikes-t); abs(spikes-t+period); ...
			       abs(spikes-t-period)]).^2)/(2 * width^2)))];
  times = [times t];
end
save('tmp');
sdf = [times; rates .* 10000/(sqrt(2 * pi) * width)];
end


function latency = sdflatency(sdf)

  pre = find(sdf(1,:) <= 300);
  prerate = mean(sdf(2,pre));
  ti = 1;
  count = 0;
  while ti < length(sdf(1,:)) & count < 5
    if sdf(2,ti) > prerate
      count = count + 1;
    else
      count = 0;
    end
    ti = ti + 1;
  end
  latency = sdf(1,ti-(count));
    
  
  