function xvars(nf, nt, resps, varargin)
%xvars(nframes, ntrials, resps, varargin)
% 'sd' s  sets the S.D. of the temporal convolution, in frames
% 'mode' 1|2  1 = convolve after noise* resp
%             2 = convolve noise, then * resp


sd = 3;
j = 1;
mode = 2;
while j < length(varargin)
    if strncmpi(varargin{j},'mode',2)
        j = j+1;
        mode = varargin{j};
    elseif strncmpi(varargin{j},'sd',2)
        j = j+1;
        sd = varargin{j};
    end
    j = j+1;
end

w=30;
k = Gauss(sd,-w:w);
r = rand(nt,nf+2*w);

if mode == 1 %convole after multiply
s = ceil(rand(nt,nf+2*w).*length(resps));
iresp = resps(s) + r.* sqrt(resps(s));
sresp = conv2(1,k,iresp,'valid');
iresp = iresp(:,w:end-w-1);
else %convolve noise first, then multiply
sr = conv2(1,k,r,'valid');
s = ceil(rand(nt,nf).*length(resps));
r= r(:,w:end-w-1);
iresp = resps(s) + r.* sqrt(resps(s));
sresp = resps(s) + sr.* sqrt(resps(s));
end
vi = var(sum(iresp'));
vs = var(sum(sresp'));
fprintf('Vars %.3f/%.3f = %.3f\n',vi,vs,vi/vs);
plot(sum(iresp'),sum(sresp'),'o');
refline(1);
corrcoef(r(:),s(:))