function fitline = plotpsych(data, mean, sd, varargin)

color = 'r';
shown = 0;
j = 1;
while(j < nargin-2)
  if(strncmpi(varargin{j},'color',5))
    j = j+1;
    color = varargin{j};
  elseif strncmpi(varargin{j},'shown',4)
    shown = 1;
  end
  j = j+1;
end

if length(data) < 2
    return;
end

h = errorbar([data.x],[data.p],sqrt([data.p] .* (1 - [data.p])./[data.n]),'o');
set(h,'color',color);
set(h,'MarkerFaceColor',color,'linewidth',2);
hold on;
step = (max([data.x]) - min([data.x]))/100;
if sum([data.x] < 0) == 0  %% xvals all positive, plot fit down to 0
x = 0:step:max([data.x]);
else
x = min([data.x]):step:max([data.x]);
end
y = 0.5 + erf((x - mean)/(sd * sqrt(2)))/2;
fitline = plot(x,y,'color',color,'linewidth',2);

if(shown)
for j = 1:length(data)
  text(data(j).x+2*step,data(j).p,sprintf('%d',data(j).n),'color',color);
end
end