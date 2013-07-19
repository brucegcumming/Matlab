function espk = ExptSpikeListAll(DATA, eid, stimes)
    t(1) = DATA.Expts{eid}.Trials(1).Start(1)-DATA.state.preperiod;
    t(2) = DATA.Expts{eid}.Trials(end).End(end)+DATA.state.postperiod;
    espk = find(stimes > t(1) & stimes < t(2));