function [value, it] = GetCheck(tag, varargin)

if strcmp(tag,'ShowSpikes') && nargin == 2
    DATA = varargin{1};
    oF = findobj(allchild(0),'flat','Tag',DATA.tag.options);
    if isempty(oF)
        value = DATA.state.showspikes;
    else
        it = findobj(oF,'Tag',tag);
        if ~isempty(it)
            value = get(it(1),'value');
        else
            value = 0;
        end
    end
    return;
end

if nargin == 2 & isfigure(varargin{1})
    it = findobj(varargin{1},'Tag',tag);
else    
    it = findobj('Tag',tag);
end
if ~isempty(it)
    value = get(it(1),'value');
else
    value = 0;
end

