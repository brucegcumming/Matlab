%Make Perl script sort out trl->mma filename, create file
%Add Gui.
%Plot Corr,AC cross sections.
%
%
%
function AC = Plot_AC_data(DPP, plottype)

  myca = [1 0 1; 0 1 0; 0.6 0.5 0; 1 0 0;  0 0 1; 0 0 0.5 ; 0.5 0 0.5; 0 0.5 0];
idx = findstr(DPP.title,'dO');
docol = 1 + (idx -1)/3;
idx = findstr(DPP.title,'ce');
cecol = 1 + (idx -1)/3;
idx = findstr(DPP.title,'SC');
sccol = 1 + (idx -1)/3;
idx = findstr(DPP.title,'st');
stcol = 1 + (idx -1)/3;
idx = findstr(DPP.title,'me');
eyecol = 1 + (idx -1)/3;


dumean = mean(DPP.data(sccol,:));
stimtypes = DPP.data(stcol,:);
eyevals = DPP.data(eyecol,:);
stim_dovals = DPP.data(docol,:);
stim_corrvals = DPP.data(cecol,:);
if(length(stim_corrvals) == 0)
  stim_corrvals = ones(size(stim_dovals));
end

for j = 1:length(DPP.strings)
  ucvals(j) = length(findstr(DPP.strings{j},'+uc'));
end


dovals = unique(stim_dovals);
cevals = [1 -1 0];

idx = find(stimtypes == DPP.mainstim & ucvals == 0 & eyevals == 0);

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
hold off;
j = 0;
strings = {};

for x = dovals;
Xrow = [];
Yrow = [];
Zrow = [];
j = j+1;
  strings = [strings {sprintf('%.2f',x)}];
  for y = cevals;
    idx = find(stim_corrvals == y & stim_dovals == x & stimtypes == ...
	       DPP.mainstim & ucvals == 0 & eyevals == 0);
    if(length(idx) > 0)
      rate = mean(DPP.data(3,idx));
      Xrow = [Xrow x];
      Yrow = [Yrow y];
      Zrow = [Zrow rate];
    end
%    fprintf('%f %f %f\n',x,y,rate);    
  end
      Xvals = [Xvals; Xrow];
    Yvals = [Yvals; Yrow];
    Zvals = [Zvals; Zrow];
end

idx = find(stimtypes == DPP.mainstim & ucvals == 1 & eyevals == 0);
if(length(idx) > 0)
      urate = mean(DPP.data(3,idx));
      plot(mean(stim_dovals(idx)),urate,'ko');
      hold on;
end

plot(Xvals(:,1),Zvals(:,1),'bo:','MarkerFaceColor','blue');
hold on;
if size(Xvals,2) > 1
  plot(Xvals(:,2),Zvals(:,2),'ro:','MarkerFaceColor','red');
end
    
title(sprintf('%s',DPP.name));


AC.disprange = [min(Xvals(:,1)) max(Xvals(:,1))];
