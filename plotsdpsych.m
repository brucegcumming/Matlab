function [plots, errs, fits] = plotsdpsych(Data,varargin)

colors = {[0 0 1], [0 1 0], [1 0 0], [1 0 1], [0 1 1], [1 1 0]};

if max(size(Data)) == 1
  plotonepsych(Data,varargin);
elseif max(size(Data)) > 1
  nplot = 1;
  for j = 1:size(Data,1)
    for k = 1:size(Data,2)
      if ~isempty(Data(j,k).resps)
	h(nplot) = plotonepsych(Data(j,k),'Color',colors{nplot});
	labels{nplot} = sprintf('%s %s',Data(j,k).label,Data(j,k).sublabel);
	nplot = nplot+1;
      end
    end
  end
  legend(h, labels,2);
end

    
  
function [ploth, errs, fits] = plotonepsych(Data,varargin)

color = 'r';
nvar = nargin - 1;
j = 1;
while j < length(varargin) & ~isempty(varargin{j})
  str = varargin{j};
  if strmatch(str,'Color')
    color = varargin{j+1};
  end
  j = j+2;
end

for j = 1:length(Data)
  h = plot(Data(j).resps(:,1),100 * Data(j).resps(:,3),'o','Color', ...
	   color);
  ploth = h;
  hold on;
  h = errorbar(Data(j).resps(:,1),100 * Data(j).resps(:,3), ...
	       100 * Data(j).resps(:,4),'o');
  set(h,'Color',color);
  errs{j} = h;
  step = (max(Data(j).resps(:,1)) -min(Data(j).resps(:,1)))/100;
  x = min(Data(j).resps(:,1)):step:max(Data(j).resps(:,1));
  sd = Data(j).fitsd * sqrt(2);
  h = plot(x,50 + 50 * erf(((x-Data(j).fitmean) ./ sd)),'Color',color);
  fits{j} = h;
end
