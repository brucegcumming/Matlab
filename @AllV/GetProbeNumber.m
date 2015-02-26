function [p, id] = GetProbeNumber(DATA, chspk, n, varargin)
% [p, id] = GetProbeNumber(DATA, chspk, id) 
%Returns probe number and index in chspk for the nth probe
if n == 3
    if length(chspk) > 2
        p = chspk(end);
        id = length(chspk);
    else
        id = 0;
        if max(chspk)  == DATA.nprobes
            p = min(chspk)-1;
        else
            p = max(chspk)+1;
        end
    end
elseif n == 4
    ispk = find(chspk == DATA.probe(1));
    if length(chspk) > 3
        if ispk ==2
            id = 3;
        elseif ispk == 3
            id = 2;
        else
            id =2;
        end
        p = chspk(id);
    end
elseif n == 1
    p = chspk(1);
    id= n;
elseif n == 2
    p = chspk(2);
    id = 2; 
end
