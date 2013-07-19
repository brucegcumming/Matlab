function res = CalcCSDPeaks(csd, varargin)
[a, maxt] = max(var(csd));
[a, maxp] = max(csd(:,maxt));
[a, minp] = min(csd(:,maxt));

    for j = 1:21
        t(j) = maxt+j-11;
    if maxp < minp
    zc = interp1(csd(maxp:minp,t(j)),maxp:minp,0);
    else
    zc = interp1(csd(minp:maxp,t(j)),minp:maxp,0);
    end
    res.zc(j) = zc;
    end

res.t = t;