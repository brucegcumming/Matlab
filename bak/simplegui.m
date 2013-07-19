function simplegui(varargin)

% Bareboes gui window with
%
% Set defaults before reading varargin
%
% name needs to be name of function - used to set callbacks etc.
name = 'simplegui';
%
%
%If strings is empty, no list is shown. Otherwise a listbox is included.
%
strings = { 'ruf915', 'ruf922'};
strings = {};

tag = [name 'Panel']; %need to change these.
init = 0;

toplevel = findobj('Tag',tag);
if ~isempty(toplevel)
  if strncmpi(varargin{1},'store',5)
    set(topelevel,'UserData',varargin{2});
    DATA = varargin{2};
  else
    DATA = get(toplevel,'UserData');
  end
end


if nargin
    if strncmpi(varargin{1},'update',5)
        update;
    elseif strncmpi(varargin{1},'close',5)
        CloseTag(tag);
    else
        j = 1;
        init = 1;
        while(j < nargin)
            if(strncmpi(varargin{j},'name',3))
                j = j+1;
                name = varargin{j};
            end
            j = j+1;
        end
    end
else
    init = 1;
end

if init & isempty(findobj('Tag',tag))
  DATA.name = name;
  InitInterface(DATA, fstrings);
end
  
function InitInterface(DATA, fstrings)

    scrsz = get(0,'Screensize');
    wsc = scrsz(3) /1280;  %scale everything as is 1280 wide
    size(1) = 380 * wsc;
    size(2) = 200 * wsc;
    listtag = [DATA.tag 'List'];
    HSPACE = 3;
    VSPACE = 2;

    cw = 9;
    ch = 18*wsc + VSPACE;
    nlines = 4;
    cntrl_box = figure('Position', [10 scrsz(4)-220*wsc 300*wsc 200*wsc],...
        'NumberTitle', 'off', 'Tag',tag,'Name',name);
    
    if( ~isempty(strings))
        lst = uicontrol(gcf, 'Style','listbox','String', strings,...
            'Callback', [name '(''setentry'')'],'Tag',listtag,...
            'Position',[10 10 size(1) size(2)-(10 + ch*nlines)]);
        
    end

    
    bp(1) = HSPACE; bp(3) = 25; bp(2) = size(2)-ch; bp(4) = 22;
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', [name '(''next'')'],...
        'String', '>>', 'Position', bp);
    bp(1) = bp(1) + bp(3) + HSPACE;
    uicontrol(gcf,'Style', 'pushbutton', 'Callback', [name '(''next'')'],...
        'String', '>>', 'Position', bp);
    
 %New row
    bp(1) = HSPACE;
    bp(2) = bp(2) - ch;
    bp(3) = 6 * cw * wsc;
    uicontrol(gcf,'Style', 'checkbox','String', 'Rewrite', 'Tag', 'ReWrite', 'Callback', [name '(''update'')'],...
        'Position', bp);
    
    bp(2) = bp(2) - ch;
    bp(3) = 6 * cw * wsc;
    uicontrol(gcf,'Style', 'text','String','Plot','Position', bp);
    bp(1) = bp(1)+bp(3)+HSPACE;
    uicontrol(gcf,'style','pop','string','CP|CEDT|Both|DT|DTCE|CP-DT|DT_CEDT|All|Extra|Spike|Seq|Psych|sdf', ...
        'Callback', [name '(''setplot'')'], 'Tag','plottype',...
        'position',bp,'value',1);
    
    bp(2) = bp(2) - ch;
    bp(1) = HSPACE;
    uicontrol(gcf,'Style', 'text','String','Text','Position', bp);
  bp(1) = bp(1)+bp(3)+HSPACE;
  uicontrol(gcf,'Style', 'edit','String',sprintf('%.2f',0.0),'Position', bp,'Tag','PopSD','Callback', ...
	    'setselect(''Update'')','Backgroundcolor',[1 1 1]);
        
    
  hm = uimenu(gcf,'Label','File');
  uimenu(hm,'Label','Close','Callback',[name '(''close'')']);
    set(gcf,'Menubar','none');
    
  set(cntrl_box,'UserData',DATA);
  
function update(varargin)
  
function CloseTag(tag)
it = findobj('Tag',tag);
if ~isempty(it)
    close(it);
end

  
  