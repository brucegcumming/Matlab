function DATA = ApplyLayout(DATA, varargin)
%DATA = ApplyLayout(DATA, varargin)
%Reads a layout file and applies to windows
%with tags labelled in DATA.tag
%uses DATA.layoutfile, or file named with
%ApplyLayout(DATA, 'layout', file
%       if file is not full pathname, then looks in
%        DATA.prefsdir (if exists)

usegui = 0;
layoutfile = [];
if isfield(DATA,'layoutfile')
    if ~exist(DATA.layoutfile)
        if isempty(fileparts(DATA.layoutfile))
           DATA.layoutfile = [DATA.gui.prefsdir '/' DATA.layoutfile '.layout.mat'];
        end
    end
    layoutfile = DATA.layoutfile;
end
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'choose',5)
        usegui = 1;
    elseif strncmpi(varargin{j},'layout',6)
        j = j+1;
        layoutfile = varargin{j};
    elseif j == 1 && ischar(varargin{j})
        layoutfile = varargin{j};
    end
    j = j+1;
end

if usegui
    if isempty(layoutfile) && isfield(DATA,'prefsdir')
        [outname, pathname] = uigetfile(DATA.prefsdir);
    else
        [outname, pathname] = uigetfile(layoutfile);
    end
    if ischar(outname)
        layoutfile = [pathname outname];
    else
        return;
    end
end

if ~exist(layoutfile)
        mycprintf('errors','Cant read Layout %s\n',layoutfile);
        return;
end

DATA.layoutfile = layoutfile;
try
    load(DATA.layoutfile);
catch
    cprintf('errors','Error loading %s\n',DATA.layoutfile);
    return;
end
    setappdata(DATA.toplevel,'Figpos',Figpos);
    f = fields(Figpos);
    for j = 1:length(f)
        it = findobj('type','figure','tag',f{j});
        if length(it) == 1
                set(it,'Position',Figpos.(f{j}));
        end
    end
    if exist('GUI','var')
        DATA.gui = CopyFields(DATA.gui, GUI, {'FontSize' 'FontName'});
    end
    if exist('Showvals','var')
        f = fields(Showvals);
        for j = 1:length(f)
            DATA.show.(f{j}) = Showvals.(f{j});
        end
    end
    if exist('Datavals','var')
        f = fields(Datavals);
        for j = 1:length(f)
            DATA.(f{j}) = Datavals.(f{j});
        end
    end
