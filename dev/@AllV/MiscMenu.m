function MiscMenu(a, b, type)    DATA = GetDataFromFig(a);    onoff = {'off' 'on'};    if sum(strcmp(type,{'savelayout' 'savedefaultlayout'}))        f = fields(DATA.tag);        for j = 1:length(f);            it = findobj('type','figure','Tag', DATA.tag.(f{j}));            if length(it) == 1                Figpos.(DATA.tag.(f{j})) = get(it,'Position');            end        end        if isempty(DATA.layoutfile) || strcmp(type,'savedefaultlayout')             DATA.layoutfile = DATA.defaultlayout;        end        [outname, pathname] = uiputfile(DATA.layoutfile);        if outname            DATA.layoutfile = [pathname outname];            save(DATA.layoutfile,'Figpos');            fprintf('Layout saved to %s\n',DATA.layoutfile);        end    elseif strcmp(type,'preferences')%        DATA = AllV.mysetappdata(DATA,'gui',DATA.gui);%        AppDataPopup(DATA, 'gui');        AllV.PreferencesPopup(DATA);    elseif strcmp(type,'loadlayout')        [afile, pathname] = uigetfile(DATA.layoutfile);        if ischar(afile)        DATA.layoutfile = [pathname afile];        AllV.ApplyLayout(DATA);        end    elseif sum(strcmp(type,{'saveconfig' 'savedefaultconfig'}))        if strcmp(type,'savedefaultconfig') || ~isfield(DATA,'configfile')            DATA.configfile = DATA.defaultconfig;        end        [outname, pathname] = uiputfile(DATA.configfile);        if outname            DATA.configfile = [pathname outname];            savefields = {'quickcutmode' 'probeswitchmode' 'auto' 'ptsz' 'plot' 'plotspk' 'gui'};            DATA.plotspk.probes = [];            SaveConfig(DATA, DATA.configfile,savefields,'verbose');        end    elseif strcmp(type,'scaledensity')        DATA.plot.scaledensity = ~DATA.plot.scaledensity;        set(a,'checked',onoff{1+DATA.plot.scaledensity});        AllV.ReplotPCs(DATA,[]);        set(DATA.toplevel,'UserData',DATA);    elseif strcmp(type,'loadconfig')        [afile, pathname] = uigetfile(DATA.configfile);        if ischar(afile)        DATA.configfile = [pathname afile];        DATA = ReadConfig(DATA, DATA.configfile);        end    elseif strcmp(type,'tofront')            FiguresToFront(DATA.tag);    elseif strcmp(type,'plotcelllist')        AllV.PlotCellList(DATA, 'showfig','reload');    endfigure(DATA.toplevel);    