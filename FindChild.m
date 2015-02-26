function it = FindChild(F, varargin)
%calls find obj for a window and any paired windows

it = findobj(F,varargin{:});
figs = findobj('type','figure');
for j = 1:length(figs)
    a = findobj(figs(j),varargin{:});
    it = [it a];
end
