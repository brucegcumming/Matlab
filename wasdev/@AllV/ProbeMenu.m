function ProbeMenu(a,b, fcn)    onoff = {'off' 'on'};[DATA, F] = GetDataFromFig(a);if strmatch(fcn,{'usecluster' 'reclassify' 'autocut' 'simple' 'reapply'})    DATA.probeswitchmode = fcn;    SetMenuCheck(a,'exclusive');    set(DATA.toplevel,'UserData',DATA);elseif strcmp(fcn,'selecttxt')    s = inputdlg('New Probe #','Change Probe',1,{num2str(AllV.ProbeNumber(DATA)+1)});    if ~isempty(s)    p = str2num(s{1});    args = {};    if DATA.fullvswitchmode.summary        args = {args{:} 'summary'};    end    AllV.ChangeProbe(DATA,[],p);    endelseif strcmp(fcn,'select')    AllV.ProbeSelector(DATA);elseif ismember(fcn,[1 2 3])    [C, DATA.Evec, DATA.pcs, DATA.dipvals, DATA.chspk] = AllV.CalcPCs(DATA,AllVoltages,fcn-1);    DATA = AllV.ReplotPCs(DATA,[]);    set(DATA.toplevel,'UserData',DATA);elseif fcn == 4    DATA.dvdt = ~DATA.dvdt;    set(a,'Checked',onoff{DATA.dvdt+1});end        