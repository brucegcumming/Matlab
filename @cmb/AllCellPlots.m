function AllCellPlots(a,b, op)
DATA = GetDataFromFig(a);
onoff = {'off' 'on'};
AllCellRes = getappdata(DATA.toplevel,'AllCellRes');
ename = Expt2Name(DATA.Expts{DATA.currentexpt(1)});
P =  getappdata(DATA.toplevel,'AllCellPlot');
if isnumeric(op)
P = cmb.NextPlot(P,AllCellRes,op);
elseif strcmp(op,'xcorr')
F = gcf;
oldname = get(F,'name');
AllExpts = getappdata(DATA.toplevel,'AllExpt');
set(F,'Name','Calculating Xcorrs');
drawnow;
xc = ExptListCorrs(AllExpts);
Exptcmb.PlotXcorr(xc,1);
set(F,'Name',oldname);
drawnow;
elseif strcmp(op,'save')
outname = [DATA.Expt.Header.fileprefix '.Cellres.' ename '.mat'];
save(outname,'AllCellRes');
fprintf('Saved Cell Results to %s\n',outname');
elseif sum(strcmp(op,{'sortbyvar' 'sortbyprobe'}))
P.(op) = ~P.(op);
set(a,'checked',onoff{1+P.(op)});

elseif strcmp(op,'next')
P = cmb.NextPlot(P,AllCellRes,0);
elseif strcmp(op,'prev')
P = cmb.NextPlot(P,AllCellRes,-1);
end
setappdata(DATA.toplevel,'AllCellPlot',P);

