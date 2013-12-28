function ChangeProbe(a,b,p)    DATA = GetDataFromFig(a);    if strcmp(p,'next')        if DATA.probelist(DATA.probe(1)) < DATA.allnprobes            p = DATA.probelist(DATA.probe(1))+1;        else            p = 0;        end    elseif strcmp(p,'prev')        if DATA.probelist(DATA.probe(1)) > 2            p = DATA.probelist(DATA.probe(1))-1;        else            p = 0;        end    elseif sum(strcmp(p,{'save' 'quicksave'}))        if strcmp(p,'save')            AllV.PCCluster(a,b,25); %save spikes and cluster        else            AllV.PCCluster(a,b,12); %save spikes and cluster        end        if DATA.probe(1) < DATA.allnprobes            p = DATA.probelist(DATA.probe(1))+1;        else            p = 0;        end    end    set(DATA.toplevel,'name',sprintf('Loading probe %d',p));    drawnow;    if ~isfield(DATA,'spkrate')        DATA.spkrate = 100;    end    if p > 0        if strncmp(DATA.DataType,'GridData',8)            if DATA.loadfromspikes                newname = DATA.name;            else                newname = regexprep(DATA.fullvname,'\.p[0-9]*FullV',sprintf('.p%dFullV',p));            end            if DATA.setspkrate > 0            AllV.AllVPcs(newname,'tchan',p,DATA.probeswitchmode,'spkrate',DATA.setspkrate);            else            AllV.AllVPcs(newname,'tchan',p,DATA.probeswitchmode);            end        else            if DATA.setspkrate > 0                AllV.AllVPcs(DATA.toplevel,'tchannew',p,DATA.probeswitchmode,'spkrate',DATA.setspkrate);            else                AllV.AllVPcs(DATA.toplevel,'tchannew',p,DATA.probeswitchmode);            end        end    end    set(DATA.toplevel,'name',get(DATA.toplevel,'Tag'));    if DATA.quickcutmode.plotspikes        DATA = get(DATA.toplevel,'UserData');    end    figure(DATA.toplevel);    