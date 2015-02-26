function C = UpdateWithNewCut(C, DATA)
%Inserts DATA.Newcut into Clusters UpdateWithNewCut(C, DATA)
%for plotting. Careful that appdata is not updated after calling this


if isfield(DATA,'NewCut') && isfield(DATA.NewCut,'probe');
    New = DATA.NewCut;
        if New.saved >= 0 && New.probe > 0 && ...
                C.exptid == DATA.NewCut.exptid && C.probe == DATA.NewCut.probe
            C.shape = New.shape;
            C.xyr = New.xyr;
            C.crit = New.crit;
        end
    end
end