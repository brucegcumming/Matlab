%Make Perl script sort out trl->mma filename, create file
%Add Gui.
%Plot Corr,AC cross sections.
%
%
%
function Plot_DPP_data(DPP, plottype)

  myca = [1 0 1; 0 1 0; 0.6 0.5 0; 1 0 0;  0 0 1; 0 0 0.5 ; 0.5 0 0.5; 0 0.5 0; 1 1 0; 1 0 1; 0 1 1];
idx = findstr(DPP.title,'dp');
dpcol = 1 + (idx -1)/3;
idx = findstr(DPP.title,'dq');
dqcol = 1 + (idx -1)/3;
idx = findstr(DPP.title,'SC');
sccol = 1 + (idx -1)/3;
idx = findstr(DPP.title,'st');
stcol = 1 + (idx -1)/3;
idx = findstr(DPP.title,'me');
eyecol = 1 + (idx -1)/3;
idx = findstr(DPP.title,'sf');
sfcol = 1 + (idx -1)/3;
idx = findstr(DPP.title,'f2');
f2col = 1 + (idx -1)/3;

dumean = mean(DPP.data(sccol,:));
stimtypes = DPP.data(stcol,:);
eyevals = DPP.data(eyecol,:);
stim_dpvals = DPP.data(dpcol,:);
stim_dqvals = DPP.data(dqcol,:);
stim_sfvals = DPP.data(sfcol,:);
for j = 1:length(DPP.strings)
  ucvals(j) = length(findstr(DPP.strings{j},'+uc'));
end


dpvals = unique(stim_dpvals);
dqvals = unique(stim_dqvals);

idx = find(stimtypes == 10 & ucvals == 0 & eyevals == 0);
sfs = [mean(DPP.data(sfcol,idx)) mean(DPP.data(f2col,idx))];

Xvals = [];
Yvals = [];
Zvals = [];
bxphases = [];
xphases = [];
yphases = [];
byphases = [];
axphases = [];
ayphases = [];
axbphases = [];
aybphases = [];
arate = [];
arateb = [];
crate = [];
crateb = [];
disps = [];
bdisps = [];
adisps = [];
abdisps = [];
j = 0;
strings = {};

for x = dpvals;
Xrow = [];
Yrow = [];
Zrow = [];
j = j+1;
  for y = dqvals;
    idx = find(stim_dqvals == y & stim_dpvals == x & stimtypes == ...
	       10 & ucvals == 0 & eyevals == 0);
    if(length(idx) > 0)
      rate = mean(DPP.data(3,idx));
      Xrow = [Xrow x];
      Yrow = [Yrow y];
      Zrow = [Zrow rate];
    end
%    fprintf('%f %f %f\n',x,y,rate);    
  
  
  end
  if length(Yrow) > 1
  Xrow = [Xrow x];
  Yrow = [Yrow (Yrow(1)+ 2 * pi)];
  Zrow = [Zrow Zrow(1)];
  else
    Yrow
  end
    disp = x/(pi * 2 * sfs(1));
    dispb = (x + 2 * pi)/(pi * 2 * sfs(1));
    iy = x * sfs(2)/sfs(1);
    iyb = (dispb * pi * 2 * sfs(2)); 
    adisp = (x - pi)/(pi *  2 * sfs(1));
    adispb = (x + pi)/(pi *  2 * sfs(1));
    ay = (adisp * pi * 2 * sfs(2)) + pi; 
    ayb = (adispb * pi * 2 * sfs(2)) + pi; 
  if(ayb > max(Yrow))
    ayb = ayb - (2 * pi);
  end
  if(ay > max(dqvals))
    ay = ay - (2 * pi);
  end
  if(ay < min(dqvals))
    ay = ay + (2 * pi);
  end
  if(iyb > max(Yrow))
    iyb = iyb - (2 * pi);
  end

  if(length(Yrow) < 2 | length(Zrow) < 2)
    Yrow
    Zrow
  else

    	guess(1) = 10;
	guess(2) = 10;
	guess(3) = 0;
  Fit.answer = nlinfit(Yrow,Zrow,'sine',guess);
  Fit.dp = x;
  Fits{j} = Fit;

  if get(findobj('Tag','ShowFit'),'value')
    iz = sine(Fits{j}.answer,iy);
   izb = sine(Fits{j}.answer,iyb);
   az = sine(Fits{j}.answer,ay);   
   azb = sine(Fits{j}.answer,ayb);
  else
   iz = interp1(Yrow,Zrow,iy);
   izb = interp1(Yrow,Zrow,iyb);
   az = interp1(Yrow,Zrow,ay);
   azb = interp1(Yrow,Zrow,ayb);
  end
   xphases = [xphases x];
   yphases = [yphases iy];
   bxphases = [bxphases x];
   byphases = [byphases iyb];
   axphases = [axphases x];
   ayphases = [ayphases ay];
   axbphases = [axbphases x];
   aybphases = [aybphases ayb];
   crate = [crate iz];
   crateb = [crateb izb];
   arate = [arate az];
   arateb = [arateb azb];
   disps = [disps disp];
   bdisps = [bdisps dispb];
   adisps = [adisps adisp];
   abdisps = [abdisps adispb];
  end
  

  if (plottype == 2)
  end
  
    Xvals = [Xvals; Xrow];
    Yvals = [Yvals; Yrow];
    Zvals = [Zvals; Zrow];
end
    Xvals = [Xvals; Xvals(1,:) + 2 * pi];
    Yvals = [Yvals; Yvals(1,:)];
    Zvals = [Zvals; Zvals(1,:)];
nplots = j;
    
xXrow = [];
xYrow = [];
xZrow = [];
for x = dpvals;
    
  idx = find(stimtypes == 3 & ucvals == 0 & eyevals == 0 & stim_sfvals == ...
	   sfs(1) & stim_dpvals == x);
    if(length(idx) > 0)
      rate = mean(DPP.data(3,idx));
      xXrow = [xXrow; x];
      xYrow = [xYrow; -(pi + 0.5)];
      xZrow = [xZrow; rate];
    end
end

yXrow = [];
yYrow = [];
yZrow = [];

for y = dqvals;
    if(y <dpvals(1))
        findy = y + 2 * pi
    else
        findy = y;
    end
    
  idx = find(stimtypes == 3 & ucvals == 0 & eyevals == 0 & stim_sfvals == ...
	   sfs(2) & abs(stim_dpvals - findy) < 0.1);
    if(length(idx) > 0)
      rate = mean(DPP.data(3,idx));
      yYrow = [yYrow y];
      yXrow = [yXrow min(min(Xvals))-0.5];
      yZrow = [yZrow rate];
    end
end
      yXrow = [yXrow min(min(Xvals)) - 0.5];
      yZrow = [yZrow rate];
      yYrow = [yYrow dqvals(1)+2*pi];
      yXrow = [yXrow min(min(Xvals))- 0.5];
      yZrow = [0 yZrow];
      yYrow = [-(pi+0.5) yYrow];
      
showfit = get(findobj('Tag','ShowFit'),'value');      
if(showfit)
  lstyle = 'o';
else
  lstyle = 'o-';
end

if(plottype == 2)
    cycle = min(min(Yvals)):0.1:max(max(Yvals));
    cycle = cycle';
    strings = {};
    hold off;
    for j = 1:length(Xvals)-1
      strings = [strings {sprintf('%.2f',Xvals(j,1))}];
%        plot(Yvals(j,:), Zvals(j,:), 'MarkerFaceColor', myca(j,:), 'Color', ...
%	     myca(j,:),'LineWidth',1,'LineStyle',lstyle,'Marker','o');
        plot(Yvals(j,:), Zvals(j,:), lstyle,'MarkerFaceColor', myca(j,:), 'Color', ...
	     myca(j,:));
    	hold on;
    end

    last = length(yYrow);
    if last >2
    plot(yYrow(2:last),yZrow(2:last),'MarkerFaceColor', myca(nplots+1,:), 'Color', ...
	     myca(nplots+1,:),'LineWidth',1,'Marker','o');
    strings = [strings {sprintf('SF=%.2f',sfs(2))}];
    end
    last = length(xXrow);
    if last >2
    plot(xXrow(2:last),xZrow(2:last),'MarkerFaceColor', myca(nplots,:), 'Color', ...
	     myca(nplots,:),'LineWidth',1,'Marker','o');
    strings = [strings {sprintf('SF=%.2f',sfs(1))}];
    end
    legend(strings);
    if(showfit)
      for j = 1:length(Xvals)-1
    	thefit = sine(Fits{j}.answer,cycle);
    	plot(cycle,thefit, 'Color', myca(j,:),'LineWidth',1);
      end
    end
    xlabel(sprintf('Phase for SF=%.2f',sfs(2)));
    ylabel('spikes');
end
if(plottype == 3)
    cycle = min(min(Xvals)):0.1:max(max(Xvals));
    cycle = cycle';
    strings = {};
    hold off;
    for j = 1:length(Yvals)-1
      strings = [strings {sprintf('%.2f',Yvals(1,j))}];
        plot(Xvals(:,j), Zvals(:,j), lstyle, 'MarkerFaceColor', myca(j,:), 'Color', ...
	     myca(j,:));
    	hold on;
    end

    last = length(xXrow);
    if last >2
    plot(xXrow(2:last),xZrow(2:last),'MarkerFaceColor', myca(nplots+1,:), 'Color', ...
	     myca(nplots+1,:),'LineWidth',1,'Marker','o');
    strings = [strings {sprintf('SF=%.2f',sfs(1))}];
    end
    last = length(yYrow);
    if last >2
    plot(yYrow(2:last),yZrow(2:last),'MarkerFaceColor', myca(nplots+2,:), 'Color', ...
	     myca(nplots+2,:),'LineWidth',1,'Marker','o');
    strings = [strings {sprintf('SF=%.2f',sfs(2))}];
    end
    
    legend(strings);
    last = size(Xvals,1);
    if(showfit)
    for j = 1:length(Yvals)-1
      Fit.answer = nlinfit(Xvals(1:last-1,j),Zvals(1:last-1,j),'sine',guess);
      Fit.dp = x;
      xFits{j} = Fit;
    	thefit = sine(xFits{j}.answer,cycle);
    	plot(cycle,thefit, 'Color', myca(j,:),'LineWidth',1);
    end
    end
    xlabel(sprintf('Phase for SF=%.2f',sfs(1)));
    ylabel('spikes');
end
save('tmp');
if (plottype == 1)
    allXvals = Xvals;
    allYvals = Yvals;
    allZvals = Zvals;
    if(length(xXrow) > 1)
      allXvals = [[xXrow; xXrow(1)+2*pi] Xvals ];
      allYvals = [[xYrow(1); xYrow] Yvals ];
      allZvals = [[xZrow(1); xZrow] Zvals ];
    end
    if(length(xXrow) > 1)
      allXvals = [yXrow; allXvals];
      allYvals = [yYrow; allYvals];
      allZvals = [yZrow; allZvals];
    end

    hold off;
    if(get(findobj( 'Tag','Shading'),'value'))
      [xi, yi] = meshgrid(linspace(min(min(Xvals)),max(max(Xvals))),linspace(min(min(Yvals)),max(max(Yvals))));
      zi = Interpf(Xvals,Yvals,Zvals,xi,yi,1,0.4);
      pcolor(xi,yi,zi);
    else
      pcolor(allXvals,allYvals,allZvals);
    end

    hold on;
    h = colorbar;
    plot(xphases,yphases,'w');
    plot(bxphases,byphases,'w');
    plot(axphases,ayphases,'w:');
    plot(axbphases,aybphases,'w:');
    xlabel(sprintf('Phase for %.2f cpd',sfs(1)));
    ylabel(sprintf('Phase for %.2f cpd',sfs(2)));
elseif (plottype ==4)
  plot([disps bdisps],[crate crateb],'o');
  hold on;
  plot([adisps abdisps],[arate arateb],'or');

  dx = -1:0.05:1;
  aphases = dx * 2 * pi  * sfs(1);
  bphases = dx * 2 * pi  * sfs(2);
  acaphases = (dx * 2 * pi  * sfs(1)) + pi;
  acbphases = (dx * 2 * pi  * sfs(2)) + pi;
  xlabel('Disparity')
  ylabel('spikes')

  idx = find(aphases > max(max(Yvals)));
  while idx
    aphases(idx) = aphases(idx) - (2 * pi);
    idx = find(aphases > max(max(Yvals)));
  end
  
  idx = find(aphases < min(min(Yvals)));
  while idx
    aphases(idx) = aphases(idx) + (2 * pi);
    idx = find(aphases < min(min(Yvals)));
  end
  
  idx = find(acaphases > max(max(Yvals)));
  while idx
    acaphases(idx) = acaphases(idx) - (2 * pi);
    idx = find(acaphases > max(max(Yvals)));
  end
  
  idx = find(acaphases < min(min(Yvals)));
  while idx
    acaphases(idx) = acaphases(idx) + (2 * pi);
    idx = find(acaphases < min(min(Yvals)));
  end
  
  idx = find(acbphases > max(max(Yvals)));
  while idx
    acbphases(idx) = acbphases(idx) - (2 * pi);
    idx = find(acbphases > max(max(Yvals)));
  end
  
  idx = find(acbphases < min(min(Yvals)));
  while idx
    acbphases(idx) = acbphases(idx) + (2 * pi);
    idx = find(acbphases < min(min(Yvals)));
  end
  
  idx = find(bphases < min(min(Yvals)));
  while idx
    bphases(idx) = bphases(idx) + (2 * pi);
    idx = find(bphases < min(min(Yvals)));
  end
  
  idx = find(bphases > max(max(Yvals)));
  while idx
    bphases(idx) = bphases(idx) - (2 * pi);
    idx = find(bphases > max(max(Yvals)));
  end
  
  irates = griddata(Xvals, Yvals, Zvals, aphases, bphases);
  airates = griddata(Xvals, Yvals, Zvals, acaphases, acbphases);
  plot(dx, irates,'b');
  plot(dx, airates,'r');
  
end
%xlabel('Phase')
if(plottype ~= 4)
  title(sprintf('%s',DPP.name));
end



