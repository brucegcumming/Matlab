function args = PlotArgs(DATA, Expt,varargin)

combined = 0;
if DATA.state.plotcombined == 1
    combined = 1;
end
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'combined',6)
        combined = 1;
    end
    j = j+1;
end

if isfield(DATA.tag,'rcfiga')
    args = {'rcfiga' DATA.tag.rcfiga 'figb' DATA.tag.rcfigb 'figc' DATA.tag.rcfigc};
else
    args = {};
end
if isempty(DATA.plot.showcp)
    DATA.plot.showcp = 0;
end
if DATA.plot.showcp > 0
    psych = 1;
else
    psych = get(findobj('Tag','PlotPsych','Parent',DATA.toplevel),'value');
end
if ~isfield(Expt.Trials,'RespDir') %% no psych here
    psych = 0;
end
seq = get(findobj('Tag','PlotSeq','Parent',DATA.toplevel),'value')-1;
on = get(findobj('Tag','ShowN','Parent',DATA.toplevel),'value');

if on
    args = {args{:} 'shown'};
end
on = get(findobj('Tag','PlotSeq','Parent',DATA.toplevel),'value')-1;
if psych & ismember(seq,[4])
    args = {args{:} 'psych' 'cpseq'};
    % noplot stops PSF from being plotted, still shows seq
    %        args = {args{:} 'psychnoplot' 'cpseq'};
else
    if psych
        args = {args{:} 'psych'};
    end
    if seq == 1
        args = {args{:} 'seqt'};
    elseif seq == 3
        args = {args{:} 'seqid'};
    elseif seq == 2
        args = {args{:} 'sequence'};
    elseif seq == 12 && combined == 1
        args = {args{:} 'sdfsub' 'result.x==0'};
    elseif seq == 13 && combined == 1
        args = {args{:} 'sdf' 'collapse' 1};
    elseif seq == 14 && combined == 1
        args = {args{:} 'sdf' 'collapse' 2};
    elseif seq == 15 && combined == 1
        args = {args{:} 'plotall'};
    elseif seq == 9 && combined == 1
        args = {args{:} 'sdfall'};
    elseif seq == 7 && combined == 1
        args = {args{:} 'cpseq'};
    elseif seq == 8 && combined == 1
        args = {args{:} 'cphist'};
    elseif seq == 17 && combined == 1
        args = {args{:} 'pcolor'};
    end
end

if ismember(DATA.plot.showcp,[1 2 3 4])
    args = {args{:},  'CPtag', DATA.tag.cptag};
end
if DATA.plot.showcp == 1
    args = {args{:} 'psych'  'cpt'};
elseif DATA.plot.showcp == 2 && combined == 0
    args = {args{:} 'cp'};
elseif DATA.plot.showcp == 3 && combined == 0
    args = {args{:} 'cphist'};
elseif DATA.plot.showcp == 4 && combined == 0
    args = {args{:} 'cpt' 'emdiff' DATA.em 'emskip' DATA.plot.emskip};
end

needsdf = 0;
for j = 1:length(args)
    if ischar(args{j}) && strncmp(args{j},'sdf',3)
        needsdf = 1;
    end
end
if needsdf
    args = {args{:} 'preperiod' DATA.state.preperiod 'postperiod' DATA.state.postperiod};
end
if DATA.plot.showem
    args = {args{:} 'eyem'};
end
if DATA.plot.addhash
    args = {args{:} 'showmu'};
end

if (DATA.plot.condenseRC || combined == 0) && isempty(strmatch(DATA.plot.type ,'Subspace'))
    args = {args{:} 'condense'};
elseif Expt.Header.rc == 2 && ~strcmp(DATA.plot.type ,'Subspace')
    args = {args{:} 'condense'};
end
if DATA.plot.centerRFmeasures
    args = {args{:} 'centerRF'};
end

if DATA.state.uselfp
    args = {args{:} 'lfpt'};
end

if DATA.plot.nminrc > 0
    args = {args{:} 'rcnmin' DATA.plot.nminrc};
end
if DATA.plot.nmin > 0
    args = {args{:} 'nmin' DATA.plot.nmin};
end
if DATA.plot.sdfw > 0
    args = {args{:} 'sdfw' DATA.plot.sdfw};
else
    args = {args{:} 'fbox'};
end

if ~DATA.plot.condenseRC && Expt.Stimvals.st == 4 && strcmp(Expt.Stimvals.et,'Ol')
    args = {args{:} 'twoslice'};
end

if strfind(Expt.Header.expname,'tfXip')
    args = {args{:} 'sxcx'};
end

if DATA.plot.collapse
    args = {args{:} 'collapse' 2};
end
if DATA.plot.flip
    args = {args{:} 'reverse'};
end
if DATA.plot.acov
    args = {args{:} 'acov'};
end
if DATA.plot.plotmod
    args = {args{:} 'mod'};
end


