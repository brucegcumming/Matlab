function res  = CheckGridV(name, varargin)
%
% CheckGridV reads in FullV files for a set of probes and checks the
% results. Default is to plot a short segment superimposing all channels. 

smpls = 7562500:7563000;
probes = 1:96;
normalize = 1;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'normalize',6)
        normalize = 0;
    end
    j = j+1;
end


A = [];
for p = probes
    pname = regexprep(name,'.p[0-9]*FullV',sprintf('.p%dFullV',p));
    load(pname);
    if normalize
        if isempty(A)
            sumv = double(FullV.V)./std(double(FullV.V));
        else
            sumv = sumv + double(FullV.V)./std(double(FullV.V));
        end
    else
        if isempty(A)
            sumv = double(FullV.V);
        else
            sumv = sumv + double(FullV.V);
        end
    end
    A(p,:) = double(FullV.V(smpls));
    if normalize
        A(p,:) = A(p,:)./std(A(p,:));
    end
end
sumv = sumv./length(probes);
    
hold off;
plot(A');
res.A = A;
res.sumv = sumv;
hold on;
plot(mean(A),'k','linewidth',3);