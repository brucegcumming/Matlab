function MakeFigure(varargin)

results = [];
resps = [];
tags{1} = 'DTsep';
fig = 1;
fsz = 12; 
j = 1;
while j <= length(varargin)
    if iscell(varargin{j})
        results = varargin{j};
    end
    j = j+1;
end

if fig == 2
load dorisfstim.mat;
subplot(2,2,1);
fillpmesh(res.disps,res.dgs,res.resps{end},'plot');
subplot(2,2,2);
imagesc(res.sresps{end});

elseif fig == 1
%
% results{5} needs to be dg x dx, for dsf and dor, as in dxdor.mat
% results{4} needs to be dori x ori, for a fixed slant as in oridori.mat
% results{3} needs to be dsf x dori, for a fixed slant, and is oridsf.mat
%
    if isempty(results)
    load dorisims.mat
    end

subplot(3,1,1);
hold off;
fillpmesh(results{5}.dgs, results{5}.disps,results{5}.resps{1},'plot');
axis('square');
%shading('interp');
subplot(3,1,2);
fillpmesh(results{5}.dgs, results{5}.disps,fliplr(results{5}.sresps{1}),'plot');
%shading('interp');
axis('square');
ylabel('disparity (degrees)','fontsiz',fsz);
xlabel('disparity gradient (degrees/degree)','fontsiz',fsz);

subplot(3,1,3);

odo= PlotSims(results{4});
oresp = [0 max(odo.sig')];
sfo = PlotSims(results{3});
subplot(3,1,3);
hold off; 
oris = [0 results{4}.oris];
plot(oris .* 180/pi, oresp,'linewidth',2);
sresp = [0 max(sfo.sig' .* -1)];
hold on;
oris = [0 results{3}.oris];
plot(oris.*180/pi, sresp,'r','linewidth',2);
set(gca,'xlim',[0 90],'xtick',[0 45 90],'ytick',[]);
axis('square');
xlabel('RF orientation (degrees)','fontsiz',fsz);
ylabel('Response','fontsiz',fsz);
end

