%Make Perl script sort out trl->mma filename, create file
%Add Gui.
%Plot Corr,AC cross sections.
%
%
%
function Plot_DPSF_data(DPP, plottype, AC)

  myca = [1 0 1; 0 1 0; 0.6 0.5 0; 1 0 0;  0 0 1; 0 0 0.5 ; 0.5 0 0.5; 0 0.5 0; 1 1 0; 1 0 1; 0 1 1];
idx = findstr(DPP.title,'dp');
dpcol = 1 + (idx -1)/3;
idx = findstr(DPP.title,'SC');
sccol = 1 + (idx -1)/3;
idx = findstr(DPP.title,'st');
stcol = 1 + (idx -1)/3;
idx = findstr(DPP.title,'me');
eyecol = 1 + (idx -1)/3;
idx = findstr(DPP.title,'sf');
sfcol = 1 + (idx -1)/3;

dumean = mean(DPP.data(sccol,:));
stimtypes = DPP.data(stcol,:);
eyevals = DPP.data(eyecol,:);
stim_dpvals = DPP.data(dpcol,:);
stim_sfvals = DPP.data(sfcol,:);
for j = 1:length(DPP.strings)
  ucvals(j) = length(findstr(DPP.strings{j},'+uc'));
end


dpvals = unique(stim_dpvals);
sfvals = unique(stim_sfvals);

idx = find(stimtypes == 3 & ucvals == 0 & eyevals == 0);

Xvals = [];
Yvals = [];
Zvals = [];

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
  for y = sfvals;
    idx = find(stim_sfvals == y & stim_dpvals == x & stimtypes == ...
	       3 & ucvals == 0 & eyevals == 0);
    if(length(idx) > 0)
      rate = mean(DPP.data(3,idx));
      Xrow = [Xrow x];
      Yrow = [Yrow y];
      Zrow = [Zrow rate];
    end
%    fprintf('%f %f %f\n',x,y,rate);    
  end
  if length(Yrow) > 1 & plottype == 1
    yinc = (Yrow(length(Yrow)) - Yrow(1))/length(Yrow);
    Xrow = [Xrow x];
    Yrow = [Yrow Yrow(length(Yrow))+yinc];
    Zrow = [Zrow Zrow(length(Zrow))];
  end
	  
  if(length(Yrow) < 2 | length(Zrow) < 2)
    Yrow
    Zrow
  else

    	guess(1) = mean(Zrow);
	guess(2) = std(Zrow);
	guess(3) = 0;
  Fit.answer = nlinfit(Yrow,Zrow,'sine',guess);
  Fit.dp = x;
  Fit.sf = y;
  Fits(j) = Fit;

  end
  
    Xvals = [Xvals; Xrow];
    Yvals = [Yvals; Yrow];
    Zvals = [Zvals; Zrow];
end

    Xvals = [Xvals; Xvals(1,:) + 2 * pi];
    Yvals = [Yvals; Yvals(1,:)];
    Zvals = [Zvals; Zvals(1,:)];
nplots = j;
    

lXrow = [];
lYrow = [];
lZrow = [];
rXrow = [];
rYrow = [];
rZrow = [];

for y = sfvals;
  idxl = find(stimtypes == 3 & ucvals == 0 & eyevals == -1 & stim_sfvals ...
	     == y);
  idxr = find(stimtypes == 3 & ucvals == 0 & eyevals == 1 & stim_sfvals ...
	     == y);
    if(length(idxr) > 0)
      rate = mean(DPP.data(3,idxr));
      rYrow = [rYrow y];
      rXrow = [rXrow min(min(Xvals))-0.5];
      rZrow = [rZrow rate];
      end
    if(length(idxl) > 0)
      rate = mean(DPP.data(3,idxl));
      lYrow = [lYrow y];
      lXrow = [lXrow min(min(Xvals))-1];
      lZrow = [lZrow rate];
    end
end

nsf = length(sfvals);
if plottype == 1  %pseudocolor - add an extra row
  if length(lXrow) > 2
    lXrow = [lXrow lXrow(nsf)];
      lZrow = [lZrow lZrow(nsf)];
      lYrow = [lYrow lYrow(nsf)+yinc];
  end
  if length(rXrow) > 2
      rXrow = [rXrow rXrow(nsf)];
      rZrow = [rZrow rZrow(nsf)];
      rYrow = [rYrow lYrow(nsf)+yinc];
  end
end

      
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
    for j = 1:length(Xvals)-1
      strings = [strings {sprintf('%.2f',Xvals(j,1))}];
%        plot(Yvals(j,:), Zvals(j,:), 'MarkerFaceColor', myca(j,:), 'Color', ...
%	     myca(j,:),'LineWidth',1,'LineStyle',lstyle,'Marker','o');
        plot(Yvals(j,:), Zvals(j,:), lstyle,'MarkerFaceColor', myca(j,:), 'Color', ...
	     myca(j,:));
    	hold on;
    end

    last = length(rYrow);
    if last >2
    plot(rYrow(2:last),rZrow(2:last),'MarkerFaceColor', myca(nplots+1,:), 'Color', ...
	     myca(nplots+1,:),'LineWidth',1,'Marker','o');
    plot(lYrow(2:last),lZrow(2:last),'MarkerFaceColor', myca(nplots+1,:), 'Color', ...
	     myca(nplots+1,:),'LineWidth',1,'Marker','o');
   strings = [strings {'Right'} {'Left'}];
    end
    legend(strings);
    if(showfit)
      for j = 1:length(Xvals)-1
    	thefit = sine(Fits(j).answer,cycle);
    	plot(cycle,thefit, 'Color', myca(j,:),'LineWidth',1);
      end
    end
    xlabel('SF');
    ylabel('spikes');
end
if(plottype == 3)
    cycle = min(min(Xvals)):0.1:max(max(Xvals));
    cycle = cycle';
    strings = {};
    for j = 1:size(Yvals,2)
      strings = [strings {sprintf('%.2f',Yvals(1,j))}];
        plot(Xvals(:,j), Zvals(:,j), lstyle, 'MarkerFaceColor', myca(j,:), 'Color', ...
	     myca(j,:));
    	hold on;
    end

    legend(strings);
    last = size(Xvals,1);
    if(showfit)
    for j = 1:length(Yvals)-1
      Fit.answer = nlinfit(Xvals(1:last-1,j),Zvals(1:last-1,j),'sine',guess);
      Fit.dp = x;
      xFits(j) = Fit;
    	thefit = sine(xFits(j).answer,cycle);
    	plot(cycle,thefit, 'Color', myca(j,:),'LineWidth',1);
    end
    end
    xlabel(sprintf('Phase Disparity'));
    ylabel('spikes');
end
save('tmp');
if (plottype == 1) %psuedocolor
    allXvals = Xvals;
    allYvals = Yvals;
    allZvals = Zvals;
    if(length(rYrow) > 1)
      allXvals = [rXrow; allXvals];
      allYvals = [rYrow; allYvals];
      allZvals = [rZrow; allZvals];
    end
    if(length(lYrow) > 1)
      allXvals = [lXrow; allXvals];
      allYvals = [lYrow; allYvals];
      allZvals = [lZrow; allZvals];
    end
save('tmp');
    if(get(findobj( 'Tag','Shading'),'value'))
      [xi, yi] = meshgrid(linspace(min(min(Xvals)),max(max(Xvals))),linspace(min(min(Yvals)),max(max(Yvals))));
      zi = Interpf(Xvals,Yvals,Zvals,xi,yi,1,0.4);
      pcolor(xi,yi,zi);
    else
      pcolor(allXvals,allYvals,allZvals);
    end

    hold on;
    h = colorbar;
    xlabel(sprintf('Phase'));
    ylabel(sprintf('SF'));
elseif (plottype == 4)
  last = size(Xvals,1);
  for j = 1:size(Yvals,2)
      Fit.answer = nlinfit(Xvals(1:last-1,j),Zvals(1:last-1,j),'sine',guess);
      Fit.dp = x;
      Fit.sf = sfvals(j);
      xFits(j) = Fit;
  end
  dx = AC.disprange(1):0.01:AC.disprange(2);

  
  for j = 1:nsf
    save('tmp');
    fitresps(:,j) = sine(xFits(j).answer,dx' .* (pi * 2 * sfvals(j)));
    acfitresps(:,j) = sine(xFits(j).answer,(dx' .* (pi * 2 * sfvals(j)) +pi));
  end
  fitresp = mean(fitresps,2);
    plot(dx,fitresp);
    hold on;
    plot(dx,mean(acfitresps,2),'r');
  xlabel('Disparity')
  ylabel('spikes')
elseif (plottype == 5)
  last = size(Xvals,1);
  for j = 1:size(Yvals,2)
      Fit.answer = nlinfit(Xvals(1:last-1,j),Zvals(1:last-1,j),'sine',guess);
      Fit.dp = x;
      Fit.sf = Yvals(1,j);
      xFits(j) = Fit;
  end
  save('tmp');
  sfs = [xFits.sf];
  fits = [xFits.answer];
  subplot(1,2,1);
  plot(sfs,fits(1,:),'ko-'); %baseline
  hold on;
  plot(sfs,abs(fits(2,:)),'bo-'); %amplitude
  plot(sfs,lZrow,'ro-');
  plot(sfs,rZrow,'go-');
  legend('Fit mean','Fit Amp','Left','Right');
  xlabel('SF')
  ylabel('Resp')
  subplot(1,2,2);
  plot(lZrow,abs(fits(2,:)),'ro');
  hold on;
  plot(rZrow,abs(fits(2,:)),'go');
  plot((rZrow + lZrow),fits(1,:),'bo');
end
%xlabel('Phase')
title(sprintf('%s',DPP.name));
save('tmp');


