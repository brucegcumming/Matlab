function [fign, new] = GetFigure(tag,varargin)

%fign =Getfigure(tag,...)   is a convenience routine that finds a figure
%whose 'Tag' property matches tag. Returns figure handle.
%If no figure is found, a new
%figure is created, thus GetFigure is guaranteed to return a valid
%figure handle
%GetFigure(tag,'front',...   brings the figure to the front.
%any additional arguments are passed on to both findobj() and figure()
%

doclear = 0;
j = 1;
fk = 1;
nk = 1;
varnew = {'Name',tag};
varfind = {};
if strcmp(computer,'PCWIN')
    forcefront = 1;
else
    forcefront = 0;
end

while j < nargin
    if strncmpi(varargin{j},'clear',4)
        doclear = 1;
    elseif strncmpi(varargin{j},'front',4)
        forcefront = 1;
    elseif strncmpi(varargin{j},'noforce',4)
        forcefront = 0;
    elseif strncmpi(varargin{j},'title',4)
        j = j+1;
        varnew{nk} = 'Name';
        varnew{nk+1} = varargin{j};
        nk = nk+2;
    else
        varnew{nk} = varargin{j};
        varfind{fk} = varargin{j};
        fk = fk+1;
        nk = nk+1;
    end
    j = j+1;
end

new = 0;
if isempty(tag)
    figure(gcf);
    return;
end
if isnumeric(tag)
    fign = tag;
else
    fign = findobj('Tag',tag,'Type','Figure',varfind{:});
end
if isempty(fign) 
    new = 1;
  fign = figure('Tag',tag,varnew{:});
elseif fign == 0
  fign = figure;
else
    if forcefront
        figure(fign);
    else
        set(0,'CurrentFigure',fign(1));
    end
end

if doclear
    delete(allchild(gca));
end
