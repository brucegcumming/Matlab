function C = CheckForMean(DATA,C)            if ismember(C.space(1),[3 4]) | (C.space(1) ==6 && C.space(2) == 4) %need mean        if ~isfield(C,'MeanSpike')            if isempty(DATA.clid) && isfield(C,'r') && isfield(C,'sign')                DATA.clid = find(C.r .* C.sign > C.crit(1) .* C.sign);                DATA.nid = setdiff(1:DATA.nevents,DATA.clid)';            end            C.MeanSpike = AllV.PlotMeanSpike(DATA,'recalc');        end    end    if length(C.mahal) > 3 & ~isfield(C,'gmdprime')        C.gmdprime = C.mahal(4);    end