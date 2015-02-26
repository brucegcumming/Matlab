function [Expt, res] = PlotCombined(DATA, Expt, varargin)
res = [];
xargs = {};
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'hold',4)
        xargs = {xargs{:} varargin{j}};
    end
    j = j+1;
end
if ~isfield(Expt.Header,'probe')
    Expt.Header.probe = DATA.probe;
end
if isfield(Expt,'ExcludeCluster')
    for j = 1:length(Expt.ExcludeCluster{1})
        id = find([Expt.Trials.Trial] == Expt.ExcludeCluster{1}(j));
        if id
            Expt.Trials(id).Trial = -Expt.Trials(id).Trial;
        end
    end
end
if ismember(DATA.plot.showcp,[5 6]) & isfield(Expt.Trials,'RespDir')
    psfargs = {};
    on = get(findobj('Tag','ShowN','Parent',DATA.toplevel),'value');
    if on
        psfargs = {psfargs{:} 'shown'};
    end
    cmb.SetFigure(DATA.tag.psych, DATA);
    hold off;
    if DATA.plot.showcp == 6
        ExptPsych(Expt, 'smooth',10,psfargs{:});
    else
        [a,b] = ExptPsych(Expt,psfargs{:});
        if isfield(b,'seedrpts')
            nbad = sum(b.seedrpts ==1);
            nseed = length(b.seedrpts);
            fraction = nbad./nseed;
            if fraction >0.2
                title(sprintf('Only %d/%d seeds repeated',nbad,nseed));
            end
        end
    end
    return;
end
cmb.SetFigure(DATA.tag.dataplot,DATA);
if sum(strcmp('Dc',{Expt.Stimvals.et Expt.Stimvals.e2})) && isfield(Expt.Trials,'Dc')
    Expt.Header.rc = 2;
elseif isfield(Expt.Stimvals,'Dc') && Expt.Stimvals.Dc > 0 && Expt.Stimvals.Dc < 1
    Expt.Header.rc = 2;
end
if Expt.Header.rc == 1
    cmb.SetFigure(DATA.tag.rcfiga,DATA);
    cmb.SetFigure(DATA.tag.rcfigb,DATA);
    cmb.SetFigure(DATA.tag.rcfigc,DATA);
end
args = cmb.PlotArgs(DATA,Expt,'combined');
args = {xargs{:} args{:}};

if isfield(DATA,'plotargs')
    args = {args{:} DATA.plotargs{:}};
end
res = PlotExpt(Expt,args{:});
if DATA.state.plotseq == 4 & isfield(res,'fig')
    set(res(1).fig, 'WindowButtonDownFcn',@cmb.FigButtonPressed);
    set(res(1).fig, 'WindowButtonUpFcn',@cmb.FigButtonReleased);
    dat = get(res(1).fig,'UserData');
    dat.parentfigtag = DATA.tag.top;
    set(res(1).fig,'UserData',dat);
elseif DATA.state.plotseq == 10  %show acloop plot
    PlotRC(res,'acloop');
elseif DATA.state.plotseq == 11  %show acloop plot
    PlotRC(res,'acresp');
end
%? SHow autocorrelation if its been calculated? 
t = get(get(gca,'title'),'String');
spikelist = cmb.WhichClusters(DATA.toplevel);
edepth = GetEval(Expt,'ed');
%    title([t 'Cl' sprintf(' %d',cmb.WhichClusters(DATA.toplevel))]);
title([t sprintf(' P%d',DATA.probe) sprintf(' Cl %d',spikelist) sprintf('ed %.2f',edepth)]);
if DATA.plot.showem
    GetFigure(DATA.tag.emplot);
    Expt = LoadEmData(Expt);
    PlotExptEM(Expt);
end
if DATA.state.uselfp
    GetFigure('LFP');
    hold off;
    CalcLFPPulse(Expt,DATA.AllData,'plot');
    if size(Expt.Trials(1).LFP,2) > 5
        if DATA.plot.lfpplot == 1
            PlotMLFP(Expt,'stack',0);
        elseif DATA.plot.lfpplot ==2
            PlotMLFP(Expt,'image');
        end
    end
end

%    save(DATA.outname,'Expt'); %use save to do the saving....


