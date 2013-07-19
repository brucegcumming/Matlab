function dcontrast(type,varargin)
%
%dcontrast(type,...)
%simulates BEM for interocular contrast differences;
%
%dcontrast('rls') does random noise (1D)
%dcontrast('grating') does gratings (also in 1D)
%
% ...,'nreps',N) set number of dot patterns to average over
%
% ...,'norm', [R R B]) sets gain of normalization for each eye, and binoc
%                      default = [0 0 0], no normalization for contrast.

%
%
% with differnt normalization gains in the two eyes, and a binocular
% normalization, response at low contrast can exceed high contrast!
% binocular normalization has a bigger effect on disparity modulation for
% RDS than for sines...
% norm [1 1 100] also shows this.  I guess product of gains > 1. Need to
% watch out for this...
nreps = 200;
ngain = [1 1 1 1];
ocu = [1 1];
j = 1;
while j < nargin
    if strncmpi(varargin{j},'norm',3)
        j = j+1;
        ngain = varargin{j};
    elseif strncmpi(varargin{j},'nreps',3)
        j = j+1;
        nreps = varargin{j};
    elseif strncmpi(varargin{j},'ocularity',3)
        j = j+1;
        ocu = varargin{j};
    end
    j = j+1;
end

GetFigure('Delta Contrast');
if strncmpi(type,'rls',3)
    [d,r] = rls(nreps,'lc',1,'norm',ngain);
    hold off;
    plot(d,r,'k');
    [d,r] = rls(nreps,'lc',0.2,'norm',ngain);
    hold on;
    plot(d,r,'r');
    [d,r] = rls(nreps,'rc',0.2,'norm',ngain);
    plot(d,r,'g');
    [d,r] = rls(nreps,'rc',0.2,'lc',0.2,'norm',ngain);
    plot(d,r,'b');
elseif strncmpi(type,'grating',3)
    np = 1;
    norm = [1 1 1];
    [p,d,r,A] = grating([ocu(1) ocu(2)],norm);
    pwr(np) = p(1);
    means(np) = p(2);
    labels{np} = sprintf('(%.2f/%.2f + %.2f/%.2f)/%.2f)',A.contrast(1),A.gains(1),A.contrast(2),A.gains(2),A.gains(3));
    np = np+1;
    hold off;
    plot(d,r,'k');
    
    norm = [1 0.2 * ngain(2)  ngain(3)];
    [p,d,r,A] = grating([ocu(1) ocu(2) .* 0.2], norm);
    pwr(np) = p(1);
    means(np) = p(2);
    labels{np} = sprintf('(%.2f/%.2f + %.2f/%.2f)/%.2f)',A.contrast(1),A.gains(1),A.contrast(2),A.gains(2),A.gains(3));
    np = np+1;
    hold on;
    plot(d,r,'r');
    
    norm = [0.2 * ngain(1) 1 ngain(3)];
    [p,d,r,A] = grating([0.2 * ocu(1) ocu(2)],norm);
    pwr(np) = p(1);
    means(np) = p(2);
    plot(d,r,'g');
    labels{np} = sprintf('(%.2f/%.2f + %.2f/%.2f)/%.2f)',A.contrast(1),A.gains(1),A.contrast(2),A.gains(2),A.gains(3));
    np = np+1;
    
    norm = [0.2 * ngain(1) 0.2 * ngain(2) ngain(4)];
    [p,d,r,A] = grating([0.2 * ocu(1) 0.2 * ocu(2)],norm);
    means(np) = p(2);
    pwr(np) = p(1);
    labels{np} = sprintf('(%.2f/%.2f + %.2f/%.2f)/%.2f)',A.contrast(1),A.gains(1),A.contrast(2),A.gains(2),A.gains(3));
    np = np+1;
    plot(d,r,'b');
    legend(labels);
%  GetFigure('GratResps');
%  plot(pwr,means,'o-');
end

