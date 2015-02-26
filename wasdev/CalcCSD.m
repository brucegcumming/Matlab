function csd = CalcCSD(R, varargin)
%csd = CalcCSD(R, ...) calculates second derivative along dimension 1
%
%csd = CalcCSD(R, 'flip') dimension 2;
%
%csd = CalcCSD(R, 'smooth', K) smooths with 2-D kernel K (sd in time, probes).
%CalcCSD(R, 'smooth') applies default smoothing (4 time samples, 2 probes);
% When used with 
%'flip' the kernal is transposed inside here...
dim = 1;

j = 1;
sk = [];
while j <= length(varargin)
    if strncmpi(varargin{j}, 'flip',4)
        dim = 2;
    elseif strncmpi(varargin{j}, 'smooth',4)
        if length(varargin) > j && isnumeric(varargin{j+1})
            j = j+1;
            sd = varargin{j};
        else %default smoothing is 2 probes and 4 timesamples
            sd = [4 2];
        end
        [X,Y] = meshgrid([-(sd(1)*2):(sd(1)*2)],[-(sd(2)*2):(sd(2)*2)]);
        sd(find(sd == 0)) = 0.1;
        sd = 1./sd; 
        sk = exp(-(X.^2./2.*sd(1).^2 + (Y.^2)./2.*sd(2).^2));
    end
    j = j+1;
end

if isfield(R,'Trials') %given and exp structure
    Expt = R;
    R = cat(3,Expt.Trials.LFP);
    R(isnan(R)) = 0;
    R = mean(R,3);
    if diff(size(R)) < 0
        R = R';
    end        
end
if dim == 2
    sk = sk';
end
csd = diff(R,2,dim);
if ~isempty(sk)
    csd = conv2(csd,sk,'same');
end