function DATA = FinishNewData(DATA,varargin)

    DATA.energy = [];
    DATA.spkvar = [];
    AllVoltages = AllV.mygetappdata(DATA,'AllVoltages');
    if isfield(DATA,'trueprobe')
        vprobe = DATA.trueprobe;
    else
        vprobe = DATA.probe(1);
    end

    DATA.energy(1,:) = sqrt(squeeze(sum(diff(AllVoltages(DATA.probe(1),:,:),1,2).^2,2)));
    xspk = DATA.chspk(DATA.chspk ~= vprobe);
    for j = 1:length(xspk)
        DATA.energy(j+1,:) = sqrt(squeeze(sum(diff(AllVoltages(xspk(j),:,:),1,2).^2,2)));
        DATA.eprobe(j+1) = xspk(j);
    end

    DATA.spkvar(DATA.probe(1),:) = squeeze(std(AllVoltages(DATA.probe(1),:,:),[],2));
    DATA.nevents = size(AllVoltages,3);
