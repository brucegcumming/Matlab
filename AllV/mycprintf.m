function mycprintf(varargin)
%wrapper for cprintf that keeps currentfigure
f = get(0,'currentfigure');
cprintf(varargin{:});
if ~isempty(f) && f > 0
    if isfigure(f)
        set(0,'currentfigure',f);
    else
        cprintf('red','!!!!Figure %.0f is no longer a figure!!!\n',f);
    end
end