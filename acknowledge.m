function acknowledge(str, tag,varargin)
%acknowledge(msg, windowlabel)  is just a wrapper for msgbox
%acknowledge(msg, windowlabel,'beep') makes a beep too 
%acknowledge(msg, windowlabel,'print') writes message to console too
%if windowlabel is a figure handle, then uses Name of figure as label
%    and set the ParentFigure appdata
% acknowledge(msg, F, 'label' windowlabel) forces the banner label
show = 0;
okstr = 'OK';
parentfig = [];
if nargin < 2
    tag = 'acknowledge';
elseif isfigure(tag)
    parentfig = tag;
    tag = get(parentfig,'name');
end
j = 1;

while j <= length(varargin)
    if strncmpi(varargin{j},'beep',4)
        beep;
    elseif strncmpi(varargin{j},'label',4)
        j = j+1;
        tag = varargin{j};
    elseif strncmpi(varargin{j},'print',4)
        show = 1;
    end
    j = j+1;
end

if show
    fprintf('%s\n',str);
end
CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle='replace';
h = msgbox(str,tag,CreateStruct);
if ~isempty(parentfig)
    setappdata(h,'ParentFigure',parentfig)
end