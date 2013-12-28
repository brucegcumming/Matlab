function [C, Evec, pcs, dip, chspk, errs, details] = CalcPCs(DATA, AllVoltages, nprobepc,  varargin)errs = {};details = [];uid = DATA.uid;Evec = [];if nargout > 3    calcfit = 1;else    calcfit = 0;endj = 1;while j <= length(varargin)    if isstruct(varargin{j}) && isfield(varargin{j},'Evec')        Evec = varargin{j};        fprintf('Using Previous EigenVectors on%s\n',sprintf(' %d',Evec.chspk));    elseif strncmp(varargin{j},'nofit',5)        calcfit = 0;    end    j = j+1;endif length(nprobepc) > 1 || nprobepc(1) >= 0    if length(nprobepc) > 1        chspk = DATA.probe(1) + nprobepc;    else        chspk = DATA.probe(1)-nprobepc:1:DATA.probe(1)+nprobepc;    end    if DATA.csd == 2        chspk = chspk-2;        chspk = chspk(chspk > 0 & chspk < size(AllVoltages,1)-1);        if isempty(chspk)            chspk = DATA.probe(1);        end    elseif DATA.csd        chspk = chspk-1;        chspk = chspk(chspk > 0 & chspk < size(AllVoltages,1));    else        chspk = chspk(chspk > 0 & chspk <= size(AllVoltages,1));    endelseif nprobepc <0     chspk = DATA.chspk;endif isfield(Evec,'chspk')    chspk = Evec.chspk;else    chspk = DATA.chspk; %with different array types, can't guess this hereendevspace = 0;if DATA.csd     if DATA.csd == 2        csd = diff(AllVoltages,2,1);        evspace = 1;    else        csd = diff(AllVoltages,1,1);        evspace = 2;    end    chspk = chspk -1;    chspk = chspk(chspk > 0 & chspk <= size(csd,1));    TV = csd(chspk(1),:,uid);    for j = 2:length(chspk)        TV = cat(2,TV,csd(chspk(j),:,uid));    endelseif DATA.dvdt == 2    TV = AllVoltages(chspk(1),:,uid);    TV = cat(2,TV,diff(AllVoltages(chspk(1),:,uid),1,2));    for j = 2:length(chspk)        TV = cat(2,TV,AllVoltages(chspk(j),:,uid));        TV = cat(2,TV,diff(AllVoltages(chspk(j),:,uid),1,2));    end    evspace = 3;elseif DATA.dvdt == 3    TV = diff(AllVoltages(chspk(1),:,uid),2,2);    for j = 2:length(chspk)        TV = cat(2,TV,diff(AllVoltages(chspk(j),:,uid),2,2));    end    evspace = 4;elseif DATA.dvdt    TV = diff(AllVoltages(chspk(1),:,uid),1,2);    for j = 2:length(chspk)        TV = cat(2,TV,diff(AllVoltages(chspk(j),:,uid),1,2));    end    evspace = 5;else    evspace = 6;    TV = AllVoltages(chspk(1),:,uid);    for j = 2:length(chspk)        TV = cat(2,TV,AllVoltages(chspk(j),:,uid));    endendTV = squeeze(TV)';if isempty(Evec) || size(TV,2) ~= size(Evec.Evec,1)     C = cov(TV);    [pc, E] = eig(C);    [a,b] = max(diag(E));    Evec.Eval = diag(E);    if b == 1        fprintf('Max Eigernvector First\n');        errs = {errs{:} 'Max Eigernvector First\n'};    else        pc = fliplr(pc); %put largest first;        Evec.Eval = Evec.Eval(end:-1:1);    end%arbitrary sign convention just so that its consistent when re-applying        for j = 1:size(pc,1)        if sum(pc(:,j)) < 0            pc(:,j) = -pc(:,j);        end    end    Evec.Evec = pc;    Evec.calctime = now;    Evec.chspk = chspk;    Evec.npts = length(uid);    Evec.space = evspace;    Evec.npts = length(uid);    Evec.forced = 0;else    Evec.forced = 1;    Evec.newspace = evspace;    C = [];endpcs = TV*Evec.Evec;if calcfit    if DATA.usebmi    for j = 1:10%        dip(j) = HartigansDipTest(sort(pcs(:,j)));        dip(j) = AllV.BimodalCoeff(pcs(:,j),1.5);    end    p = DATA.pcplots;    for j =1:length(p)        [as(j),bs(j)] = AllV.BestAngle(pcs(:,p(j,1)),pcs(:,p(j,2)),1);    end    dip = bs;    [P, dip(length(p)+1)] = AllV.GMfit(pcs(:,1:4),2,1);    else    [P, dip(1), details] = AllV.GMfit(pcs(:,1:4),2,1);    endelse    dip = NaN;end