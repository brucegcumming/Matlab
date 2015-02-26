function DATA = DoLayout(a,b, fcn)

if nargin == 2
    fcn = b;
end
DATA = GetDataFromFig(a);
if strcmp(fcn,'choose')
    DATA = ApplyLayout(DATA, 'choose');
elseif strcmp(fcn,'load')
    if ~exist(DATA.layoutfile)
        name = DATA.layoutfile;
        d = dir([DATA.gui.prefsdir '*.layout.mat']);
        id = find(CellToMat(strfind({d.name},name)));
        if ~isempty(id)
            if length(id) > 1
                fprintf('possible layout files:\n');
                for j = 1:length(id)
                    fprintf('%s\n',d(id(j)).name);
                end
                [a,b] = max([d(id).datenum]);

                fprintf('Using %s\n',d(id(b)).name);
                DATA.layoutfile = [DATA.gui.prefsdir '/' d(id(b)).name];
            else
                DATA.layoutfile = [DATA.gui.prefsdir '/' d.name];
            end
        end
    end
    if exist(DATA.layoutfile)
        DATA = ApplyLayout(DATA,DATA.layoutfile);
    end
elseif strcmp(fcn,'loadlast')
    name = strrep(DATA.defaultlayout,'.layout','last.layout');
    DATA = ApplyLayout(DATA, 'layout', name);
elseif strcmp(fcn,'setfont')
    font = uisetfont();
    DATA.gui = copyfields(DATA.gui,font,{'FontSize' 'FontName' 'FonwWeight'});
    SetUIRecursive(DATA.toplevel,font);
    set(DATA.toplevel,'UserData',DATA);
elseif strcmp(fcn,'savelast')
    name = strrep(DATA.defaultlayout,'.layout','last.layout');
    DATA = SaveLayout(DATA, 'layout', name,'nochoose');
elseif strcmp(fcn,'savedefault')
    DATA = SaveLayout(DATA, 'layout', DATA.defaultlayout);
elseif strcmp(fcn,'save')
    DATA = SaveLayout(DATA);
end


