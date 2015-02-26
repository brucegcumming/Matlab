function SetProbeList(DATA)
it = findobj(DATA.toplevel,'Tag','ProbeId');
if isfield(DATA,'ArrayConfig') && isfield(DATA.ArrayConfig,'id') && length(DATA.probenames) == length(DATA.ArrayConfig.id)
[a,b] = sort(DATA.ArrayConfig.id);
set(it,'string',DATA.probenames(b),'value',DATA.probe);
setappdata(it,'probelist',b);
else
set(it,'string',DATA.probenames,'value',DATA.probe);
end


