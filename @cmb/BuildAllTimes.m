function DATA = BuildAllTimes(DATA, varargin)
%Make sure spkcahce is up to data for  for all expts

expts = 1:length(DATA.Expts);
probes = 1:length(DATA.probelist);

j= 1;
while j <= length(varargin)
    if strcmp(varargin{j},'expts')
        j = j+1;
        expts = varargin{j};
    end
    j = j+1;
end

DATA = cmb.SpkCache(DATA,max(expts),max(probes),'init');
OLD.currentexpt = DATA.currentexpt;
OLD = CopyFields(OLD, DATA, 'probe');
for e = expts(:)'
    %DATA.currentexpt = e;
    for p = 1:length(probes)
        if cmb.SpkCache(DATA,e,p,'check') == 0
            DATA.probe = p;
            DATA = cmb.LoadSpikes(DATA,e);
            DATA.spklist = ExptSpikeListAll(DATA, e, DATA.AllData.Spikes.times);
            DATA = CalcClusterVars(DATA,  DATA.spklist,'expt',e);
            DATA = cmb.SetCluster(DATA, e, p);
            [DATA, DATA.spklist] = SetExptSpikes(DATA,DATA.currentexpt(1),0,'useexpt');
        end
    end
end
DATA = CopyFields(DATA,OLD);
%set(DATA.toplevel,'UserData',DATA);
