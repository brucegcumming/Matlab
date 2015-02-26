function FullV = bin2FullV(bin, chan, name)

FullV.V(1,:) = double(bin(chan,:))./32000;
if nargin < 3
FullV.name = 'C:/smr/david/David';
else
    FullV.name = name;
end
smoothw = 40;

FullV.samper = 1./24000;
FullV.exptno = 1;
FullV.chspk = chan;
if 0
    ts = now;
    sm = smooth(double(FullV.V(1,:)),smoothw);
    FullV.V = FullV.V - sm;
    FullV.filtertime = mytoc(ts);
end
FullV.intscale = [2 32000];
