function args = CopyUiProperties(it, varargin)
%args = CopyUiProperties(it)
%returns {'property' get(it,'property') for a list of properties
%default list is label, tag, callback, style,
%used to create a duplicate uicontrol:
%args = CopyUiProperties(btn);
%uicontrol(F, args{:}) makes a bopy of btn in F
args = {};
otype = get(it,'type')
if strcmp(otype,'uimenu')
    p = {'label' 'tag' 'callback'};
else
    p = {'label' 'tag' 'callback' 'style'};
end
j = 1;
while j <= length(varargin)
    if iscell(varargin{j})
        p = varargin{j};
    elseif ischar(varargin{j})
        p = {p{:} varargin{j}};
    end
    j = j+1;
end

for j = 1:length(p)
    args = {args{:} p{j} get(it,p{j})};
end
   

