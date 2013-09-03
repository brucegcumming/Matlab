function DATA = SaveLayout(DATA, varargin)
%DATA = ApplyLayout(DATA, varargin)
%Reads a layout file and applies to windows
%with tags labelled in DATA.tag
%uses DATA.layoutfile, or file named with
%ApplyLayout(DATA, 'layout', file
%       if file is not full pathname, then looks in
%        DATA.prefsdir (if exists)

layoutfile = [];
if isfield(DATA,'layoutfile')
    layoutfile = DATA.layoutfile;
elseif isfield(DATA,'defaultlayout')
    layoutfile = DATA.defaultlayout;    
end

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'layout',5)
        j = j+1;
        layoutfile = varargin{j};
    elseif j == 1 && ischar(varargin{j})
        layoutfile = varargin{j};
    end
    j = j+1;
end


savefields = {'Figpos'};
if isfield(DATA,'gui')
    GUI = DATA.gui;
    savefields = {savefields{:} 'GUI'};
end

Figpos = getappdata(DATA.toplevel,'Figpos');
f = fields(DATA.tag);
   for j = 1:length(f);
       it = findobj('type','figure','Tag', DATA.tag.(f{j}));
       if length(it) == 1
           Figpos.(DATA.tag.(f{j})) = get(it,'Position');
       end
   end
   [outname, pathname] = uiputfile(layoutfile);
   if outname
       DATA.layoutfile = [pathname outname];
       save(DATA.layoutfile,savefields{:});
       fprintf('Layout saved to %s\n',DATA.layoutfile);
   end
   setappdata(DATA.toplevel,'Figpos',Figpos);
   