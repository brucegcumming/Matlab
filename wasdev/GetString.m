function f = GetString(tag, parent, callback,varargin)
%f = GetString(tag, parent, callback)
%pops up a simple window with one text box.
%sends string to named callback

lbl = 'Done';
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'label',5)
        j=j+1;
        lbl = varargin{j};
    end
    j=j+1;
end

[f,isnew] = GetFigure(tag);
if isnew
    x = get(parent,'Position');
    set(f,'position',[x(1) x(2) 300 50],'menubar','none');
    uicontrol(f, 'string',[],'style','edit','units','normalized','position',[0.05 0.55 0.9 0.45],...
        'Callback',callback);
    uicontrol(f, 'string','Done','style','pushbutton','units','normalized','position',[0.05 0.05 0.2 0.45],...
        'Callback',callback);
    uicontrol(f, 'string',lbl,'style','text','units','normalized','position',[0.25 0.05 0.7 0.45],...
        'Tag','CommentLabel');
    set(f,'UserData',parent);
end