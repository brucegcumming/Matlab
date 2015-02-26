function SxCx(Expt, varargin)
%caluclate modulation in Autocorrelatin functions to classify Simple/Compex
%
%
res = PlotExpt(Expt,'acov','noplot');
for j = 1:size(res.period,1)
for k = 1:size(res.period,2)
    ts = 1:res.period(j,k);
    Fn(j,k,1) = famp(ts,squeeze(res.acov(j,k,ts)),1/res.period(j,k));
    Fn(j,k,2) = famp(ts,squeeze(res.acov(j,k,ts)),2/res.period(j,k));
end
end
plot(Fn(:,:,2)+Fn(:,:,1),Fn(:,:,2)./Fn(:,:,1),'o');