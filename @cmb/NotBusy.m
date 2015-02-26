function NotBusy(DATA)

[a,b,c] = fileparts(DATA.datafilename);
if isfield(DATA,'toplevel') && ishandle(DATA.toplevel)
set(DATA.toplevel,'Name',[DATA.tag.top ':' b]);
end

