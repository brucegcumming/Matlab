function [zn,zp] = corlsim(factor,varargin)


cfifty = 0.5;
nr = 1;
nl = 1;

j = 1;
while j < nargin
    if strncmpi(varargin{j},'C50',3)
        j = j+1;
        cfifty = varargin{j};
    elseif strncmpi(varargin{j},'normalize',4)
        j = j+1;
        nl = varargin{j};
        j = j+1;
        nr = varargin{j};
    end
    j = j+1;
end
[L,R] = meshgrid(0:0.1:1);
SL = L;
NL = L *nl;
NR = R *nr;
L = L .* factor;

z = (L + R).^2;
zn = ((L - R).^2) ./ (cfifty + (NL.^2 + NR.^2));
GetFigure('CORLA');
pcolor(SL,R,zn);
colorbar;
title('NULL disp');
shading('flat');
ca = caxis;

GetFigure('CORLB');
zp = ((L + R).^2) ./ (cfifty + (NL.^2 + NR.^2));
pcolor(SL,R,zp);
colorbar;
title('PREF disp');
shading('flat');
cb = caxis;

caxis([min([ca cb]) max([ca cb])]);
colorbar;
GetFigure('CORLA');
caxis([min([ca cb]) max([ca cb])]);
colorbar;