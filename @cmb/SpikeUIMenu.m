function SpikeUIMenu(a,b);
DATA = GetDataFromFig(a);
tag = get(a,'Tag');
SetMenuCheck(a,'exclusive',1);

if strcmp(tag,'ns5')
    DATA.state.usensx = 1;
elseif strcmp(tag,'nev')
    DATA.state.usensx = 2;
else
    DATA.state.usensx = 0;
end
if strcmp(tag,'Trig-')
    DATA.triggersign = -1;
elseif strcmp(tag,'Trig+')
    DATA.triggersign = 1;
elseif strcmp(tag,'TrigBoth')
    DATA.triggersign = 0;
end
set(DATA.toplevel,'UserData',DATA);



