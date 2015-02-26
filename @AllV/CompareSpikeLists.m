function DATA = CompareSpikeLists(DATA, ida, idb,iida, iidb)
AllVoltages = AllV.mygetappdata(DATA, 'AllVoltages');
p = AllV.ProbeNumber(DATA);
oldid = find(DATA.oldclst == 2);
olnid = find(DATA.oldclst == 1);

tb = DATA.oldtrigtimes(oldid(iidb));
ta = DATA.t(DATA.clid(iida));
if isempty(tb)
    goodid = ones(size(ta));
else
    for j = 1:length(ida)
        [c,d] = min(abs(ta(j)-tb));
        tadiff(j) = c;
        if c > 0.0005
            goodid(j) = 0;
        else
            goodid(j) = 1;
        end
    end
    ida = ida(goodid == 0);
    iida = iida(goodid == 0);
end
goodid = [];
if isempty(ta)
    tbdiff = NaN;
else
    for j = 1:length(idb)
        [c,d] = min(abs(tb(j)-ta));
        if c > 0.0005
            goodid(j) = 0;
        else
            goodid(j) = 1;
        end
        tbdiff(j) = c;
    end
end
idb = idb(goodid == 0);
Vall = AllV.mygetappdata(DATA,'Vall');
GetFigure('Changed Spikes');
subplot(1,1,1);
if length(iida) > 1000
    Vs = squeeze(AllVoltages(p,:,DATA.clid(iida(1:1000))));
else
    Vs = squeeze(AllVoltages(p,:,DATA.clid(iida)));
end
hold off;
if ~isempty(Vs)
    h = plot(Vs,'r');
    h = h(1);
else
    h(1) = 0;
end
if length(iidb) > 10000
    xx = iidb(1:1000);
else
    xx = iidb;
end
%ocid is memmbers of AllVoltages that were spikes but are not any longer
ocid = find(ismember(round(DATA.t.*10000),round(DATA.oldtrigtimes(oldid(xx)).*10000))); 
oid = find(ismember(round(Vall.t.*40000),round(DATA.oldtrigtimes(oldid(xx)).*40000)));
if ~isempty(oid)
    vpts = bsxfun(@plus,oid,DATA.spts');
    Vs = double(reshape(Vall.V(p,vpts),size(vpts)));
    if isfield(Vall,'intscale')
        Vs = Vs.* Vall.intscale(1)./Vall.intscale(2);
    end
    Vs = Vs + diff(minmax(Vs(:)))/2';
    hold on;
    xx = plot(Vs,'k');
    plot(squeeze(AllVoltages(p,:,ocid))+diff(minmax(Vs(:)))/2','m');
    DATA.oldclusterpts{1} = [ocid(:)' iida(:)']; %includecurrent spikes not classified before
    DATA.oldclusterpts{1} = [ocid(:)']; %just events that were spikes
    h(2) = xx(1);
    mylegend(h, {sprintf('%d new',length(iida)) sprintf('%d old',length(oid))});
else
    fprintf('Missing SPike times not in New Vall\n');
end
if ~isempty(ocid)
    AllV.ReplotPCs(DATA,[]);
end
