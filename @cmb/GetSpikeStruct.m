function Spks=GetSpikeStruct(DATA)

if isempty(DATA.AllData.Spikes)
Spks = DATA.AllSpikes{DATA.probe};
else
Spks = DATA.AllData.Spikes;
end

