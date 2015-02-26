function ExportPlots(a,b, caller)
% ExportPlots(a,b, caller)
% Exports variables from a gui to the workspace


DATA = GetDataFromFig(a);
if strcmp(caller,'combine')
    LFPresult = getappdata(DATA.toplevel,'LFPresult');
    Expt = DATA.Expt;
    varnames = {'LFPresult' 'Expt'};
    labels = {'LFPresult' 'Expt'};
    items = {LFPresult  Expt};
    t = sprintf('Export from combine');
    export2wsdlg(labels, varnames, items, t);
end
