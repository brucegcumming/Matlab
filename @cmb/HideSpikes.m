function HideSpikes(a,b)
DATA = GetDataFromFig(a);
x = cmb.GetCheck('NoSpikes')
DATA.state.nospikes = ~x;
cmb.SetCheck('NoSpikes',DATA.state.nospikes);
if x
set(a,'String','Show Waveforms');
else
set(a,'String','Hide Waveforms');
end
set(DATA.toplevel,'UserData',DATA);


