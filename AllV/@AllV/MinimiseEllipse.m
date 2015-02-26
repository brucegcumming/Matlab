function [SSD, details ] = MinimiseEllipse(params, DATA, state)%[SSD, details ] = AllV.MinimiseEllipse(params, DATA, state)%find ellipze that maximizes separation of 2 Gaussian fits in 1dcx = params(1);cy = params(2);rx = params(3);ry = params(4);a = params(5);xy = DATA.xy{state.cluster};xys = xyrotate(xy(:,1)-cx,(xy(:,2)-cy) ./state.aspectratio,a);r = ((xys(:,1))./rx).^2 + ((xys(:,2))./(ry./state.aspectratio)).^2;%r = AllV.Rprime(r);details.r = r;C = DATA.cluster;C.clst(r < 1) = 2;C.clst(r>=1) = 1;[dp, fits] = AllV.Fit2Gauss(C, r, DATA);%don't just fit on traditional dprime.  If noise is far away, can%pay to find samll sd for cluster, even by putting boundary in the%middle of the clusterif isempty(fits{1}) || isempty(fits{2})    SSD = 1e10;elsedpa = (1 - fits{1}.mean)./fits{1}.sd;dpb = (fits{2}.mean-1)./fits{2}.sd;details.fits = fits;%don't just sum. Favor the smaller number, If one dp is large, don' t let improvements%in that swamp the smaller oneif state.mintype == 1    dp = -(sqrt(dpa)+sqrt(dpb)).^2;else    dp = -(dpa+dpb);endif dpa < 0 || dpb < 0    SSD = 0;else    SSD = dp; endendfitparams = getappdata(DATA.toplevel,'fitparams');fitparams = [fitparams; [params dp]];setappdata(DATA.toplevel,'fitparams',fitparams);rhist = histc(r,linspace(0, max(r)));fithists = getappdata(DATA.toplevel,'fithists');fithists = cat(2,fithists,rhist);setappdata(DATA.toplevel,'fithists',fithists);