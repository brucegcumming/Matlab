function ShowLFPPlot(DATA)
b = [];

cmb.SetFigure(DATA.tag.lfpplot, DATA);
a = [];
LFP = getappdata(DATA.toplevel,'LFPExpt');
if ~isfield(LFP,'Trials')
    return;
end
if ismember(DATA.plot.lfpplot,[4 5 6 7 8 9 10 11 12 13 14 15])
    if LFP.Header.rc
        if LFP.Stimvals.n2 > 4 & LFP.Stimvals.nt > 4
            [a,b] = PlotRevCorAny(LFP,'lfp',DATA.probe,'nmin',20,'figa',DATA.tag.rcfiga,'figb',DATA.tag.rcfigb,'figc',DATA.tag.rcfigc,'twoslice');
        else
            [a,b] = PlotRevCorAny(LFP,'lfp',DATA.probe,'nmin',20,'figa',DATA.tag.rcfiga,'figb',DATA.tag.rcfigb,'figc',DATA.tag.rcfigc);
        end
        a.Header = b.Header;
        GetFigure(DATA.tag.lfpplot);
    else
        DATA.probelist=1:24;
        a = PlotLFP(LFP);
        a.Header = LFP.Header;
        a.Header.exptype = regexprep(a.Header.expname,'.*\.','');
    end
elseif ismember(DATA.plot.lfpplot,[16])
    [a,b] = PlotMLFP(LFP,'lines','ftpwr');
    setappdata(DATA.toplevel,'LFPresult',b);
    fprintf('getappdata(%d, ''LFPresult'') to Get Plotted result\n',DATA.toplevel);
end
hold off;
probelist = DATA.probelist;
probelist=1:24;
if DATA.plot.lfpplot == 1
    if LFP.Header.rc > 0
        if isempty(b)
            [a,b] = PlotRevCorAny(LFP,'lfp',DATA.probe,'figa',DATA.tag.rcfiga,'figb',DATA.tag.rcfigb,'figc',DATA.tag.rcfigc,'nmin',DATA.plot.nminrc);
            a.Header = b.Header;
            GetFigure(DATA.tag.lfpplot);
        end
        PlotAllProbes(a);
    elseif strmatch('square.co',LFP.Header.expname)
        a = PlotMLFP(LFP,'image');
    else
        a = PlotLFP(LFP);
        a.Header = LFP.Header;
        a.Header.exptype = regexprep(a.Header.expname,'.*\.','');
        subplot(1,1,1);
        PlotAllProbes(a,'LFPTrial','probes',probelist);
    end
end
issdf = isfield(a,'sdfs');
if DATA.plot.lfpplot == 3
    PlotMLFP(LFP,'image','probes',probelist);
elseif DATA.plot.lfpplot == 2
    PlotMLFP(LFP,'stack','probes',probelist);
elseif DATA.plot.lfpplot == 4
    PlotAllProbes(a,'blank','probes',probelist);
elseif DATA.plot.lfpplot == 5
    
elseif DATA.plot.lfpplot == 6
    PlotAllProbes(a,'probes',probelist);
elseif DATA.plot.lfpplot == 7 & issdf
    PlotAllProbes(a,'onestim','xvals',1:length(a.sdfs.x(:,1)),'yvals',1:length(a.sdfs.x(1,:)),'probes',probelist);
elseif DATA.plot.lfpplot == 8
    PlotAllProbes(a,'monoc','probes',probelist);
elseif DATA.plot.lfpplot == 9
    PlotAllProbes(a,'LFPEig','probes',probelist);
elseif DATA.plot.lfpplot == 10
    PlotAllProbes(a,'stimvar','probes',probelist);
elseif DATA.plot.lfpplot == 11
    PlotAllProbes(a,'varblank','probes',probelist);
elseif DATA.plot.lfpplot == 12
    PlotAllProbes(a,'frameresp','probes',probelist);
elseif DATA.plot.lfpplot == 13
    PlotAllProbes(a,'LFPTrial','probes',probelist);
elseif DATA.plot.lfpplot == 14
    PlotAllProbes(a,'CSD','probes',probelist);
elseif DATA.plot.lfpplot == 15
    PlotAllProbes(a,'LFPTrialStart',[30 100],'probes',probelist);
end


