function sumdata(alldata, ranges, varargin)

setnmin = 10;
mindiff = 0.00001;
sdlist = unique([alldata.sd]);
side = 'L';
colors = mycolors;
j = 1;
while j < nargin-1
    if strncmpi(varargin{j},'nmin',4)
        j = j + 1;
        setnmin = varargin{j};
    elseif strncmpi(varargin{j},'side',4)
        j = j + 1;
        side = varargin{j};
    elseif strncmpi(varargin{j},'sds',3)
        j = j + 1;
        sdlist = varargin{j};
    elseif strncmpi(varargin{j},'mindiff',4)
        j = j + 1;
        mindiff = varargin{j};
    end
    j = j+1;
end


%Now combine across days.

suml = [];
sumr = [];
sumData = [];
ndat = 1;

j = 0;
for sd = sdlist;
  idx = find([alldata.sd] == sd);
  j = j+1;
  for val = unique([alldata.x])
    id = find([alldata.x] == val & [alldata.sd] == sd);
    if ~isempty(id)
      sumData(ndat).x = val;
      sumData(ndat).expno = j;
      sumData(ndat).n = sum([alldata(id).n]);
      sumData(ndat).resp = sum([alldata(id).resp]);
      sumData(ndat).p = sumData(ndat).resp/sumData(ndat).n;
      sumData(ndat).name = alldata(ndat).name;
      sumData(ndat).sd = sd;
      sumData(ndat).xo = mean([alldata(id).xo]);
      ndat = ndat+1;
    end
  end
end

hold off;
summ = [];

hold off;
nplot = 1;
h = [];
labels = {};
for sd = sdlist
  sddata = sumData(find([sumData.sd] == sd));
  if(abs(sd) > 19)
    xl = -1;
    xu = 1;
  else
    xl = -0.1;
    xu = 0.1;
  end
  nmin = setnmin;
  id = find([ranges.sd] == sd);
    if ~isempty(id)
      xl = ranges(id).min;
      xu = ranges(id).max;
      nmin = ranges(id).nmin;
    end
  if length(sddata) >1
    fit = fitpsf(sddata,'xmin',xl,'xmax',xu,'nmin',nmin,'sderr');
    [xs, xid] = sort([fit.data.x]);
    idx = find(diff(xs) < mindiff);
    for j = idx;
        fit.data(j).x = mean([fit.data(j:j+1).x])
        fit.data(j).n = sum([fit.data(j:j+1).n]);
        fit.data(j).resp = sum([fit.data(j:j+1).resp]);
        fit.data(j).p = fit.data(j).resp/fit.data(j).n;
        fit.data(j+1).n = 0;
        fit.data(j+1).x = NaN;
    end
    if ~isnan(fit.fit(1))
    h(nplot) = plotpsych(fit.data,fit.fit(1),fit.fit(2),'color',colors{nplot},'shown');
    labels{nplot} = sprintf('Sd = %.0f %d Trials %.3g',sd, ...
			    sum([fit.data.n]),fit.fit(2));
    drawnow;
    summ(nplot).delay = sd;
    summ(nplot).sd = fit.fit(2);
    summ(nplot).sdlim = fit.sdlim;
    summ(nplot).pse = fit.fit(1);
    summ(nplot).data = fit.data;
    summ(nplot).xo = mean([sddata.xo]);
    nplot = nplot+1;
    else
      fprintf('Error in %d\n',mean([sddata.sd]));
    end;
    hold on;
  end
end
  legend(h,labels,2);

