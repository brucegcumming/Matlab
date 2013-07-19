function RalfResample(ndata,ncells,nsample)

%ndata is number of samples per data point
%ncells is the number of data points,
%nsample is the number of cycls of resampling
if nargin<1, ndata=5; end
if nargin<2, ncells=1000; end
if nargin<3, nsample=100; end

 
ydata=randn(ndata,ncells);

 
ymean=mean(ydata);

 
disp(['var of means                         : ' num2str(var(ymean))]);

 
index=1:ndata;
for i=1:ncells
  for j=1:nsample
    [dummy idx]=Draw_Random(index);
    yres(j,i)=mean(ydata(idx,i));
  end
end

 
disp(['var of resampled means               : ' num2str(var(yres(:)))]);

 
% additional correction since I'm not recentering here
yres_cor=yres*sqrt(ndata/(ndata-1)); 

 
disp(['var of ind-corrected resampled means : ' num2str(mean(var(yres_cor)))]);

 
% correction due to higher variance of resamples
yres_cor_pop=yres_cor/sqrt((2*ndata-1)/(ndata-1)); 

 
disp(['var of pop-corrected resampled means : ' num2str(var(yres_cor_pop(:)))]);

 

 

 
%-------------------------------------------------------------

 
  function [idx_random indices] = Draw_Random(index)
    n=length(index);
    indices=floor(rand(1,n)*n+1);
    idx_random = index(indices);
  end

 
end

