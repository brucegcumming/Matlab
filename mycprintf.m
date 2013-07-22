function mycprintf(varargin)
%wrapper for cprintf that keeps currentfigure
f = get(0,'currentfigure');
cprintf(varargin{:});
if ~isempty(f) && f > 0
    set(0,'currentfigure',f);
end