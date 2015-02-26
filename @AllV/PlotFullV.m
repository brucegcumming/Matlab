function [Vall, id] = PlotFullV(DATA, t, varargin)    addmean = 0;    dohold = 0;    markprobes = 0;    color = 'k';    j = 1;    chspk = AllV.UseProbeList(DATA,DATA.plotspk.nfullvprobes);    while j <= length(varargin)        if strncmpi(varargin{j},'addmean',7)            addmean = 1;        elseif strncmpi(varargin{j},'color',5)            j = j+1;            color = varargin{j};        elseif strncmpi(varargin{j},'chspk',5)            j = j+1;            chspk = varargin{j};        elseif strncmpi(varargin{j},'hold',4)            dohold = 1;        elseif strncmpi(varargin{j},'markprobes',5)            markprobes = 1;        end        j = j+1;    end    Vall = AllV.mygetappdata(DATA,'Vall');    if size(Vall.V,1) > 1        voff = DATA.voffset- DATA.voffset(AllV.ProbeNumber(DATA));    else        voff = 0;    end    id = find(Vall.t > t(1)  & Vall.t < t(2));    if isempty(id)        fprintf('No data between %.2f and %.2f\n',t);    end    if isfield(Vall,'intscale') && isinteger(Vall.V)        vscale = Vall.intscale(1)./Vall.intscale(2);    else        vscale = 1;    end    if dohold        hold on;    else        hold off;    end    for p = 1:length(chspk)        V = double(Vall.V(chspk(p),id)).*vscale+voff(chspk(p));        if addmean && isfield(Vall,'meanV')            V = V + Vall.meanV(id);        end        plot(Vall.t(id),V,'color',color);        hold on;        if markprobes            text(Vall.t(id(1)), voff(chspk(p)), sprintf('%d',chspk(p)));        end    end    axis('tight');        