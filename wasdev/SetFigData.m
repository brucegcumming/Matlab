function SetFigData(DATA)

if isfield(DATA,'toplevel') && DATA.toplevel > 0
    set(DATA.toplevel,'UserData',DATA);
end