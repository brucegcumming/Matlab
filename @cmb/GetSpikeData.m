function Spks = GetSpikeData(DATA, p)
if isfield(DATA,'AllSpikes')
Spks = DATA.AllSpikes{p};
else
Spks = DATA.AllData.Spikes;
end

