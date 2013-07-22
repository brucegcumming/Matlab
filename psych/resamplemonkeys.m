clear all;
load dufleft.mat
for j = 1:length(suml)
  fit = fitpsf(suml(j).data);
  fit = resamplepsf(fit, 1000);
  suml(j).sdci = fit.sdci;
end
save('dufleft.mat','suml');

clear all;
load dufright.mat
for j = 1:length(sumr)
  fit = fitpsf(sumr(j).data);
  fit = resamplepsf(fit, 1000);
  sumr(j).sdci = fit.sdci;
end
save('dufright.mat','sumr');

clear all;
load rufleft.mat
for j = 1:length(suml)
  fit = fitpsf(suml(j).data);
  fit = resamplepsf(fit, 1000);
  suml(j).sdci = fit.sdci;
end
save('rufleft.mat','suml');

clear all;
load rufright.mat
for j = 1:length(sumr)
  fit = fitpsf(sumr(j).data);
  fit = resamplepsf(fit, 1000);
  sumr(j).sdci = fit.sdci;
end
save('rufright.mat','sumr');

