function X = GetDataStruct(DATA, f)                if DATA.interactive >= 0        X = getappdata(DATA.toplevel,f);    else        X = DATA.(f);    end