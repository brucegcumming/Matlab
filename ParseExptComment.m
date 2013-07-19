function Expt = ParseExptComment(Expt,s)
%Expt = ParseExptComment(Expt,s)
%read comment string s, and update any necessary values in Expt
%
if strmatch('cm=back',s)
    sscanf(s,'cm=back=%s,wi=%f,hi=%f');
    c = findstr(s,',');
    Expt.Stimvals.backstim = s(9:c(1)-1);
    b = sscanf(s(c(1):end),',%2s=%*f');
    a = sscanf(s(c(1):end),',%*2s=%f');
    for j = 1:length(a)
        f = ['back' b(j*2-1) b(j*2)];
        Expt.Stimvals.(f) = a(j);
    end
    if 0 %old method. Fails if the format changes...
    a = sscanf(s(c(1):end),',wi=%f,hi=%f,ce=%f,xo=%f,yo=%f,dx=%f,co=%f');
    Expt.Stimvals.backwi = a(1);
    Expt.Stimvals.backhi = a(2);
    Expt.Stimvals.backxo = a(4);
    Expt.Stimvals.backyo = a(5);
    Expt.Stimvals.backdx = a(6);
    Expt.Stimvals.backco = a(7);
    end
end