function ptsize = CheckPtSize(DATA, nspk)
    if DATA.plot.setptsize
        ptsize = DATA.plot.setptsize;
    elseif ~isfield(DATA,'ptsize')
        ptsize = 1;
    elseif nspk < 2000
        ptsize = 10;
    elseif nspk < 5000
        ptsize = 6;
    else
        ptsize = 1;
    end
