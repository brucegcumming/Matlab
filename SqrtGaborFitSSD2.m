function SSQ = SqrtGaborFitSSD2(params, xx,sqrtyy,sqrtyuncorr,maxamp,maxfreq,minoffset,maxoffset)


fity = Gabor(params, 'xv',xx);
SSQ = sum((sqrtyy - fity).^2);
if params(5) > maxoffset || params(5) < minoffset
    SSQ = SSQ.*2;
end
if params(1) > maxfreq
    SSQ = SSQ.* (params(1)/maxfreq)^2;
end
if params(4) > maxamp
    SSQ = SSQ.* (params(4)/maxamp)^2;
end