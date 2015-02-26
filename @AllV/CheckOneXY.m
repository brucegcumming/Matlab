function CheckOneXY(DATA, varargin)

F = GetFigure(DATA.tag.tmplscore,'findonly');
if double(F) > 0
    AllV.AddParameterMenu(F,@AllV.XYplot,'X');
    AllV.AddParameterMenu(F,@AllV.XYplot,'Y');
end
