function DATA = SaveCallback(DATA, a, varargin)

if ishandle(a)
    if isfield(DATA,'guistate')
        DATA.guistate.lastguihandle = a;
    else
        DATA.gui.lastguihandle = a;
    end
    SetData(DATA);
end