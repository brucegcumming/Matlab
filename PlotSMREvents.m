function PlotSMREVents(name)

ENDSTIM=3;
FRAMESIGNAL=5;
BADFIX=11;
STARTTRIAL=12;
STARTEXPT=1;
ENDEXPT=2;
STARTSTIM=6;
WURTZLATE=7;
WURTZOK=8;
WURTZPREM=9;
WURTZOKW=20;
ENDTRIAL=16;
CANCELEXPT=19;
BWISREADY=4;

STOREBIT=16;

load(name);
stimlvl = Ch36;
Ev = Ch31;

sid = find(Ev.codes(:,1) == STARTSTIM);
eid = find(Ev.codes(:,1) == ENDSTIM);
badid = find(Ch30.codes(:,1) == BADFIX);
nev = length(stimlvl.level);
onoff(1:2:nev*2) = stimlvl.level;
onoff(2:2:nev*2-1) = stimlvl.level(2:end);
onoff(nev*2) = 0;
onofftimes(1:2:nev*2) = stimlvl.times;
onofftimes(2:2:nev*2) = stimlvl.times;
estimes = stimlvl.times(stimlvl.level ==0);
bstimes = stimlvl.times(stimlvl.level ==1);

hold off;

plot(onofftimes,onoff);
hold on;
plot(Ev.times(sid),1.1,'ro');

for j = 1:length(badid)
    id = find(estimes > Ch30.times(badid(j)));
    bfdelay(j) = estimes(id(1))-Ch30.times(badid(j));
end
[a,b] = sort(bfdelay);

goodin = find(Ev.codes(:,1) == WURTZOK);
goodout = find(Ch30.codes(:,1) == WURTZOK);
for j = 1:length(goodin)
    id = find(estimes < Ev.times(goodin(j)));
    offdelay(j) = Ev.times(goodin(j))-estimes(id(end));
    okdelay(j)  = Ch30.times(goodout(j))-Ev.times(goodin(j));
end
lateid = find(Ev.codes(:,1) == WURTZLATE);
for j = 1:length(lateid)
    id = find(estimes < Ev.times(lateid(j)));
    oid = find(Ev.times(goodin) < Ev.times(lateid(j)));
    lateoffdelay(j) = Ev.times(goodin(oid(end)))-estimes(id(end));
    lateoffdelay(j) = Ev.times(lateid(j))-estimes(id(end));
    oid = find(Ev.times(goodin) > Ev.times(lateid(j)));
    lateoffdelay(j) = Ev.times(goodin(oid(1)))-Ev.times(lateid(j));
end

plot(offdelay,okdelay,'o');
