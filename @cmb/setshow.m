function cntrl_box = setshow(DATA, tag)
wsc = DATA.wsc;
SPACE = 3 * wsc(1);
VSPACE = 5 * wsc(2);
h = 220 * wsc(2);
w = 350 * wsc(1);
scrsz = get(0,'Screensize');
cw = DATA.plot.cw;
ch = DATA.plot.ch;
SpkDefs;
bp(1) = SPACE;
bp(2) = ch+VSPACE;
bp(3) = cw*9;
bp(4) = ch+VSPACE;
dat.parentfigtag = DATA.tag.top;
cntrl_box = figure('Position', [200 scrsz(4)-(h+30)*wsc(2) w*wsc(1) h*wsc(2)], 'Menubar', 'none',...
    'NumberTitle', 'off', 'Tag',tag,'Name','Showvals','UserData',dat);
Figpos = getappdata(DATA.toplevel,'Figpos');
if isfield(Figpos,tag)
    set(cntrl_box,'position',Figpos.(tag));
end
setappdata(cntrl_box,'ParentFigure',DATA.toplevel);

bp(1) = SPACE;
bp(2) = VSPACE;
bp(3) = cw*8;
bp(4) = ch*1.5;
    uicontrol(gcf,'Style', 'text','string','Add Field','Position', bp);
bp(1) = bp(1)+bp(3);    
bp(3) = cw*8;
    uicontrol(gcf,'Style', 'Edit','Position', bp,...
        'Tag','addfield','Callback',{@cmb.ShowUpdate, 'addfield'});


bp(2) = bp(2)+ch*1.4;
bp(1) = SPACE;
fn = fields(DATA.show);
for j = 1:length(fn)
    id = strmatch(fn{j},CodeNames.Codes);
    if length(id)
        str = CodeNames.Label{id};
    else
        str = fn{j};
    end
    uicontrol(gcf,'Style', 'CheckBox','String',str,'Position', bp,...
        'Tag',fn{j},'Callback',@cmb.ShowUpdate,'value',DATA.show.(fn{j}));
    bp(2) = bp(2)+ch*1.4;
end

