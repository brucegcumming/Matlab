%Make Perl script sort out trl->mma filename, create file
%Add Gui.
%Plot Corr,AC cross sections.
%
%
%


DPP = Read_DPP_data('/usr/data/dufus/299/duf299.0.2grating.DPP.trl');

dumean = mean(DPP.data(4,:));
stim_dpvals = DPP.data(16,:);
stim_dqvals = DPP.data(17,:);
dpvals = unique(stim_dpvals);
dqvals = unique(stim_dqvals);

Xvals = [];
Yvals = [];
Zvals = [];
for x = dpvals;
Xrow = [];
Yrow = [];
Zrow = [];
  for y = dqvals;
    idx = find(stim_dqvals == y & stim_dpvals == x);
    rate = mean(DPP.data(3,idx));
    Xrow = [Xrow x];
    Yrow = [Yrow y];
    Zrow = [Zrow rate];
%    fprintf('%f %f %f\n',x,y,rate);    
  end
    Xvals = [Xvals; Xrow];
    Yvals = [Yvals; Yrow];
    Zvals = [Zvals; Zrow];
end
shading('interp');
pcolor(Xvals,Yvals,Zvals);


