function ShowUpdate(a,b, varargin)

DATA = GetDataFromFig(a);
j = 1;
while j <= length(varargin)
    if strcmp(varargin{j},'addfield')
        DATA = AddField(DATA, a);
    end
    j = j+1;
end

fn = fields(DATA.show);
for j = 1:length(fn)
    DATA.show.(fn{j}) = cmb.GetShow(DATA,fn{j});
end
set(DATA.toplevel,'UserData',DATA);


function DATA = AddField(DATA, a)

F = GetFigure(a);
SpkDefs; %Makes CodeNames
cw = DATA.plot.cw;
ch = DATA.plot.ch;
str = get(a,'string');
id = find(strcmp(str,CodeNames.Codes));
if isempty(id)
    if ~isfield(DATA.show,str)
        DATA.show.(str) = 1;
    end
    xf = findobj(allchild(F),'flat','style','CheckBox');
    bp(2) = (1+length(xf)) * ch * 1.4;

    bp(1) = 5;
    bp(3) = DATA.plot.cw*length(str);
    bp(4) = ch;
    uicontrol(gcf,'Style', 'CheckBox','String',str,'Position', bp,...
        'Tag',str,'Callback',@cmb.ShowUpdate,'value',DATA.show.(str));
end
