function SetXYCluster(a,b, type, plot, val,varargin)
DATA = GetDataFromFig(a);
recalcxy = 0;
if type == 1
DATA.plot.clusterX = val;
recalcxy = 1;
elseif type == 2
DATA.plot.clusterY = val;
recalcxy =1;
elseif type == 3
DATA.plot.clusterZ = val;
elseif type == 4
DATA.plot.clusterX = val(1);
if length(val) > 1
DATA.plot.clusterY = val(2);
else
DATA.plot.clusterY = varargin{1};
end
elseif type == 5
DATA.plot.clusterX = val(1);
if length(val) > 2
DATA.plot.clusterY = val(2);
DATA.plot.clusterZ = val(3);
else
DATA.plot.clusterY = varargin{1};
DATA.plot.clusterZ = varargin{2};
end
elseif type == 6
DATA.clusterErange = [1:31] + (plot-1) * 31;
end
ispk = DATA.spklist;
if (plot == 3 | type == 3 ) & DATA.plot.clusterZ > 0
cmb.Plot3DClusters(DATA,recalcxy);
else
DATA = CalcClusterVars(DATA, DATA.Expts{DATA.currentexpt(1)}.gui.spks,'force');
hold off;
DATA = cmb.DrawXYPlot(DATA, ispk); %returned DATA has chnges to clusterX/Yrange
end
set(DATA.toplevel,'UserData',DATA);
if sum(ismember([DATA.plot.clusterY DATA.plot.clusterX DATA.plot.clusterZ],[33 34 37 38]))
cmb.PlotEig(DATA);
end


cmb.SetGui(DATA);


