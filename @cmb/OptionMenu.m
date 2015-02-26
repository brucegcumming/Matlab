function DATA = OptionMenu(a,b,fcn, varargin)

DATA = GetDataFromFig(a);
if ishandle(a)
    tag = get(a,'tag');
else
    tag = b;
end

if strcmp(fcn,'QuickConfig')
    d = dir([DATA.gui.prefsdir '/*.quickconfig']);
    id = find(CellToMat(strfind({d.name},tag)));
    if ~isempty(id)
        name = [DATA.gui.prefsdir '/' d(id(1)).name];
        DATA = ReadConfig(DATA, name);
    end
    cmb.SetGui(DATA);
end
  

