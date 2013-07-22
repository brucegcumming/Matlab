function TwoPulse(varargin)

%simulate effect of two motion pulses on detection.
x = 0:500;
t = [10 1];
showj =1;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'show',4)
        j = j+1;
        showj = varargin{j};
    elseif strncmpi(varargin{j},'tc',2)
        j = j+1;
        t = varargin{j};
    end
    j = j+1;
end

ta = exp(-x/t(1));
tb = exp(-x/t(2));
ta = ta./sum(ta);
tb = tb./sum(tb);
tab = fliplr(conv(fliplr(ta),tb));
tab = tab(500:499+length(ta));
plot(x,ta);
pulse = zeros(3,100);
pulse(1,50:62) = 1;
ipulse = pulse(3,:);
ipulse(5) = 1;
GetFigure('Kernels');
hold off; 
plot(ta);
hold on;
plot(tab,'r');
GetFigure('Convolutions');
hold off;

    delays = [10:2:50];
    vid = 1:size(pulse,2);
for j = 1:length(delays)
    pulse(3,delays(j):delays(j)+12) = 0.3;
    pulse(2,delays(j):3:delays(j)+12) = 0.9;
    stim = pulse(1,:) + pulse(3,:);
    resp = conv(stim,ta);
    s(j) = max(resp);
    ustim = conv(pulse(2,:),tb);
    ustim = ustim(vid)+pulse(1,:);
    ustim = conv(ustim,ta);
    pulse(3,delays(j):delays(j)+12) = 0;
    u(j) = max(ustim);
    if ismember(j, showj)
    plot(resp);
    hold on;
    plot(ustim,'r');
    end
    pulse(2,delays(j):delays(j)+12) = 0;
end
GetFigure('MaxResps');
hold off;
plot(s);
hold on;
plot(u,'r');
