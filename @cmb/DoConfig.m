function DATA = DoConfig(a,b, fcn)
DATA = GetDataFromFig(a);

if nargin== 2
    fcn = b;
end

if strcmp(fcn,'load')
    if ischar(b) && ~isempty(b)
        DATA.configfile = b;
    end
    if exist(DATA.configfile)
        DATA = ReadCombineConfig(DATA, DATA.configfile);
    else
        name = DATA.configfile;
        d = dir([DATA.gui.prefsdir '*.config']);
        id = find(CellToMat(strfind({d.name},name)));
        if ~isempty(id)
            if length(id) > 1
                fprintf('possible config files:\n');
                for j = 1:length(id)
                    fprintf('%s\n',d(id(j)).name);
                end
                [a,b] = max([d(id).datenum]);

                fprintf('Using %s\n',d(id(b)).name);
                DATA.configfile = [DATA.gui.prefsdir '/' d(id(b)).name];
            else
                DATA.configfile = [DATA.gui.prefsdir '/' d.name];
            end
            DATA = ReadCombineConfig(DATA, DATA.configfile);
        end
    end
elseif strcmp(fcn,'select')
    DATA = ReadCombineConfig(DATA, 'choose');
elseif strcmp(fcn,'savelast')
    name = strrep(DATA.defaultconfig,'.config','last.config');
    SaveConfig(DATA, name, {'state' 'plot' 'show' 'options'});
elseif strcmp(fcn,'savedefault')
    SaveConfig(DATA, DATA.defaultconfig, {'state' 'plot' 'show' 'options'}, 'choose');
elseif strcmp(fcn,'save')
    DATA.configfile = SaveConfig(DATA, DATA.configfile, {'state' 'plot' 'show'}, 'choose');
end
SetData(DATA);


function DATA = ReadCombineConfig(DATA, name)
%Call ReadConfig, but reset some things that are really
%state dependent, like state.online
D = DATA;
DATA = ReadConfig(DATA, name);
DATA.state.online = D.state.online;
DATA.state.nospikes = D.state.nospikes; 
