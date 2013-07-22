function ddi(varargin)
%simulate properties of DDI. 
%This looks at effect of scaling tuning curves on DDI, calculated with and
%without the sqrt transform....


maxs = [10:10:100];
mins = [0:10:90];
fano = 1;
scale = 0.5;
pmaxs = maxs.*scale;
pmins = mins.*scale;
for j = 1:length(maxs)
    for k = 1:j
        preddi(j,k) = (sqrt(maxs(j))-sqrt(mins(k)))./(sqrt(maxs(j))-sqrt(mins(k)) + 2.*fano);
        postddi(j,k) = (sqrt(pmaxs(j))-sqrt(pmins(k)))./(sqrt(pmaxs(j))-sqrt(pmins(k)) + 2.*fano);
        rawpreddi(j,k) = ((maxs(j))-(mins(k)))./(maxs(j)-mins(k) + (sqrt(fano.*maxs(j)) + sqrt(fano .*mins(j))));
        rawpostddi(j,k) = ((pmaxs(j))-(pmins(k)))./(pmaxs(j)-(pmins(k)) + (sqrt(fano.*pmaxs(j)) + sqrt(fano .*pmins(j))));
    end
end

hold off;
plot(preddi(:),postddi(:),'o');
hold on;
plot(rawpreddi(:),rawpostddi(:),'ro');
