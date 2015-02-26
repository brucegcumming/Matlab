function PlotExptsSummary(Expts, varargin)
%PlotExptsSummary(Expts) Summarizes expts in a cell array or PsychFile 

plottype = 'hbar';
if ischar(Expts)
    name = Expts;
    Expts = ReadPsychFile(name, 'useallexpts');
end
Comments = {};

j = 1;
while j <= length(varargin)
    if isstruct(varargin{j}) && isfield(varargin{j},'comment')
        Comments = varargin{j};
    end
    j = j+1;
end

if isempty(Expts)
    return;
end
x = [];
y = [];
for j = 1:length(Expts)
    if isfield(Expts{j},'Trials')
    starts(j) = ConvertTime(Expts{j},Expts{j}.Trials(1).Start(1));
    if isfield(Expts{j}.Trials,'End')
        ends(j) = ConvertTime(Expts{j},Expts{j}.Trials(end).End(end));
    else
        ends(j) = ConvertTime(Expts{j},Expts{j}.Trials(end).Start(end)) +2/(24 *60*60);
    end
    names{j} = Expt2Name(Expts{j});
    x(end+1) = starts(j);
    x(end+1) = ends(j);
    nt(j) = length(Expts{j}.Trials);
    y(end+1:end+2) = length(Expts{j}.Trials);
    stored(j) = 1; %presumably, unless this is online
    if isfield(Expts{j},'stored')
        stored(j) = Expts{j}.stored;
    end
    end
end

durs = ends-starts;
enames = unique(names);
if strcmp(plottype,'hbar')
    dots = plot([nt nt(end)],[starts ends(end)],'.');
    types = unique(names);
    labelled = zeros(size(types));
    colors = mycolors;
    datetick('y');
    for j = 1:length(starts)
        e = find(strcmp(names{j},types));
        h = patch([0 0 nt(j) nt(j)],[starts(j) ends(j) ends(j) starts(j)], colors{e},'edgecolor',colors{e}); 
        if stored(j) ~= 1
            set(h,'facecolor','none');
        end
        if isfield(Expts{j},'Comments')
            cid = find(strncmp('cm=',Expts{j}.Comments.text,3));
            xid = find(strncmp('cm=rf=',Expts{j}.Comments.text,6));
            xid = [xid find(strncmp('cm=back',Expts{j}.Comments.text,7))];
            xid = [xid find(strncmp('cm=obback',Expts{j}.Comments.text,9))];
            cid = setdiff(cid,xid);
            for  c = cid(:)'
                y = ConvertTime(Expts{j},Expts{j}.Comments.times(c));
                text(nt(j),y,Expts{j}.Comments.text{c}(4:end),'horizontalalignment','right');
            end
        end
    end
    for j = 1:length(types)
        id = find(strcmp(types{j},names));
        [a,b] = max(durs(id));
        a = id(b);
        text(0,(starts(a)+ends(a))/2,names{a},'verticalalignment','middle','fontsize',14);
    end
    xlabel('Number of Trials');
    ylabel('Time');
    fname = GetName(Expts);
    if iscellstr(fname);
        [a,b] = Counts(fname);
        [c,d] = max(a);
        fname = b{d};
    end
    title(sprintf('%s on %s',fname,datestr(ConvertTime(Expts{1},0),'dd-mmm-yyyy')));
    delete(dots);
    xl = get(gca,'xlim');
    for j = 1:length(Comments)
        text(xl(2),Comments(j).date,Comments(j).comment,'HorizontalAlignment','Right','fontsize',14);
    end
else
    plot(x,y);
end