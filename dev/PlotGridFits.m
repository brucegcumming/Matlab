function result =PlotGridFits(Expts, varargin)

result = [];
plottype = 'Image';
pvarcrit = 0.3;
ArrayConfig = [];

j = 1;
while j <= length(varargin)
if strncmpi(varargin{j},'Array',4)
    j = j+1;
    ArrayConfig = varargin{j};
elseif strncmpi(varargin{j},'pvarcrit',4)
    j = j+1;
    pvarcrit = varargin{j};
elseif strncmpi(varargin{j},'GridPos',6)
    plottype = varargin{j};
end
j = j+1;
end

[X,Y] = meshgrid(1:10,1:10);
Z(10,10) = NaN;
Z(1,1) = NaN;
Z(1,10) = NaN;
Z(10,1) = NaN;
rfpos = zeros(length(Expts),2);
if ~isempty(ArrayConfig)
    for j = 1:length(ArrayConfig.X)
        zid(j) = find(X == ArrayConfig.X(j) & Y == ArrayConfig.Y(j));
        DATA.ProbeMesh(ArrayConfig.X(j),ArrayConfig.Y(j)) = j;
    end
else
    zid = [2:9 11:90 92:99];
    DATA.ProbeMesh = Z;
    DATA.ProbeMesh(zid) = [1:length(zid)];
end
DATA.ProbeMesh = DATA.ProbeMesh';

for j = 1:length(Expts)
    if Expts{j}.fit(1).pvar > pvarcrit && Expts{j}.fit(2).pvar > pvarcrit
        rfpos(j,:) = Expts{j}.fit(1).RFpos;
        rfsz(j,1) = Expts{j}.fit(1).sd;
        rfsz(j,2) = Expts{j}.fit(2).sd;
    else
        rfpos(j,:) = NaN;
    end
end
GetFigure('GridPlot');
clf
if strcmp(plottype,'Image')
subplot(2,1,1);
Z(zid) = rfpos(:,1);
Z(DATA.ProbeMesh == 0) = NaN;
[x,y,z] = fillpmesh(X,Y,Z);
h = pcolor(x,y,z);
set(h,'ButtonDownFcn', @HitImage);
shading('flat');
colorbar;

subplot(2,1,2);
Z(zid) = rfpos(:,2);
[x,y,z] = fillpmesh(X,Y,Z);
h = pcolor(x,y,z);
set(h,'ButtonDownFcn', @HitImage);
shading('flat');
colorbar;
elseif strcmp(plottype,'GridPos')
  for j = 1:size(X,1)
      nx = 1;
      x = [];
      y = [];
      plist = [];
      for k = 1:size(X,2)
          p = DATA.ProbeMesh(j,k);
          if ~isnan(p) && p > 0 && ~isnan(rfpos(p,1))
              x(nx) = rfpos(p,1);
              y(nx) = rfpos(p,2);
              plist(nx) = p;
              nx = nx+1;
          end
      end
      plot(x,y,'r-');
      hold on;
  end
  for j = 1:size(X,2)
      nx = 1;
      x = [];
      y = [];
      plist = [];
      for k = 1:size(X,1)
          p = DATA.ProbeMesh(k,j);
          if ~isnan(p) && p > 0 && ~isnan(rfpos(p,1))
              x(nx) = rfpos(p,1);
              y(nx) = rfpos(p,2);
              plist(nx) = p;
              nx = nx+1;
          end
      end
      plot(x,y,'b-');
      hold on;
  end
for j = 1:size(X,2)
    for k = 1:size(X,1)
        p = DATA.ProbeMesh(k,j);
        if ~isnan(p) && p > 0 && ~isnan(rfpos(p,1))
            plot(rfpos(p,1),rfpos(p,2),'ro','buttondownfcn',{@HitPoint, p},'markerfacecolor','r');
        end
    end
end
  id = find(~isnan(rfpos(:,1)));
  rfmean = mean(rfpos(id,:))
  result.rfpos = rfpos(id,:);
  result.rfsz = rfsz(id,:);
end


DATA.Expts = Expts;
DATA.toplevel = gcf;
set(DATA.toplevel,'UserData',DATA);

function HitImage(a,b)

pos = get(gca,'currentpoint');
DATA = get(gcf,'UserData');
xy = floor([pos(1,1) pos(1,2)]);
p = DATA.ProbeMesh(xy(2),xy(1));
GetFigure('GridOnePlot');
hold off;
PlotResult(DATA.Expts{p}.plotres);
F = DATA.Expts{p}.fit;
fprintf('At %.0f,%.0f (P%d), %s\n',floor(pos(1,1)),floor(pos(1,2)),p,F(1).fitstr);

function HitPoint(a,b, p)

pos = get(gca,'currentpoint');
DATA = get(gcf,'UserData');
GetFigure('GridOnePlot');
hold off;
PlotResult(DATA.Expts{p}.plotres);
F = DATA.Expts{p}.fit;
fprintf('At %.0f,%.0f (P%d), %s\n',floor(pos(1,1)),floor(pos(1,2)),p,F(1).fitstr);


