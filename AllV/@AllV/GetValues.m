function [x, xid] = GetValues(DATA, name, dim)    x = [];    xid = 0;    tid = find(strcmp(name,DATA.TemplateLabels));    if isempty(tid)        if sum(strncmp(name,{'ADC', 'Cen'},3))            AllVoltages = AllV.mygetappdata(DATA,'AllVoltages');        end        if strncmp(name,'PC',2)           pc = sscanf(name,'PC%d');           x = DATA.pcs(:,pc);           xid = pc;        elseif strcmp(name,'ADC1csd')            AllVoltages = diff(AllVoltages,2,1);            p = DATA.vpts(1,1);            if p > 1                p = p-1;            elseif p > size(AllVoltgates,1)                p = size(AllVoltages,1);            end            x = squeeze(AllVoltages(p,DATA.vpts(1,2),:));        elseif strncmp(name,'ADC1dvdy',8)            AllVoltages = diff(AllVoltages,1,1);            p = DATA.vpts(1,1);            if p > 1 && strcmp(name,'ADC1dvdy1')                p = p-1;            elseif p > size(AllVoltages,1)                p = size(AllVoltages,1);            end            x = squeeze(AllVoltages(p,DATA.vpts(1,2),:));                    elseif strncmp(name,'ADC',3)            if strfind(name,'dvdt')                AllVoltages = diff(AllVoltages,[],2);            end            if strncmp(name,'ADC1',4)                x = squeeze(AllVoltages(DATA.vpts(1,1),DATA.vpts(1,2),:));                xid = DATA.vpts(1,1);            elseif strncmp(name,'ADC2',4)                x =squeeze(AllVoltages(DATA.vpts(1,3),DATA.vpts(1,4),:));                xid = DATA.vpts(1,1);            elseif regexp(name,'ADC [0-9]*:')                [a,b] = sscanf(name,'ADC %d:%d');                x =squeeze(AllVoltages(a(1),a(2),:));                xid = a(1);            end        elseif strcmp(name,'energysum')            x = sum(DATA.energy)';        elseif strcmp(name,'energy')            x = DATA.energy';        elseif regexp(name,'energy[0-9]+')            p = sscanf(name,'energy%d');            p = find(DATA.chspk == p);            if size(DATA.energy,1) >= p                x = DATA.energy(p,:)';            end        elseif strcmp(name,'spkvar')            x = DATA.spkvar';        elseif strcmp(name,'CentroidSplit')%Area under poisitive region post spike            t = find(DATA.spts ==0);            C = squeeze(AllVoltages(DATA.probe(1),t:end,:));            C(C < 0) = 0;            x = sum(C.^2)';        elseif strcmp(name,'CentroidSplit')%centroid of positiv points only                        C = squeeze(AllVoltages(DATA.probe(1),:,:));            C(C < 0) = 0;            x = C' * DATA.spts'./sum(C)';        elseif strcmp(name,'Centroid')            C = squeeze(AllVoltages(DATA.probe(1),:,:)).^2;            x = C'.^2 * DATA.spts'./sum(C.^2)';        else            if strcmp(name,'rTchan')                xid = 1;            elseif dim(1) == 3                xid = dim(2);            end            if ~isempty(xid) &&  isfield(DATA,'TemplateScores')                x = DATA.TemplateScores(:,xid);                            end        end    else        if isempty(tid)        else            if isfield(DATA,'TemplateScores')            x = DATA.TemplateScores(:,tid(1));            xid = tid(1);            else                x = DATA.pcs(:,5);            end        end    end            