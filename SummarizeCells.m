function [Expts, S] = SummarizeCells(path, varargin)

S = [];
compareplots = [];

j = 1;
while j <= length(varargin)
    if strcmpi(varargin{j},'compare')
        j = j+1;
        compareplots = varargin{j};
    end
    
    j = j+1;
end

if isstruct(path) && isfield(path,'type') %previous result
    Expts = path;
elseif isdir(path)
    d = dir([path '/*.cell*.mat']);
    nx = 0;
    
    for j = 1:length(d)
        if strfind(d(j).name,'Cells')
        else
            nx = nx+1;
            name = [path '/' d(j).name];
            Expt = LoadExpt(name);
            Expts(nx).type = Expt.Header.expname;
            if isfield(Expt,'fit') && isfield(Expt.fit,'peak');
                Expts(nx).peak = Expt.fit.peak;
            else
                Expts(nx).peak = NaN;
            end
            if isfield(Expt,'fit') && isfield(Expt.fit,'sd');
                Expts(nx).sd = Expt.fit.sd;
            else
                Expts(nx).sd = NaN;
            end
            if isfield(Expt.Header,'cellnumber')
                Expts(nx).cell = Expt.Header.cellnumber;
            end
            if isfield(Expt.Header,'probe')
                Expts(nx).probe= Expt.Header.probe;
            end
        end
    end
end


types = unique({Expts.type});
for j = 1:length(types)
    id = strmatch(types{j},{Expts.type},'exact');
    S(j).type = types{j};
    for k = 1:length(id)
        if isempty(Expts(id(k)).peak)
            Expts(id(k)).peak = NaN;
        end
    S(j).probe(k) = Expts(id(k)).probe ;
    S(j).cell(k) = Expts(id(k)).cell;
    S(j).peak(k) = Expts(id(k)).peak;
    S(j).sd(k) = Expts(id(k)).sd;
    end
end



if compareplots
    CompareExpts(S(compareplots(1)),S(compareplots(2)));
end


function CompareExpts(Ea,Eb)

hold off;
[id, idb] = ismember(Ea.cell,Eb.cell);
for j = 1:length(id)
    if id(j)
        or = [Ea.peak(j) Eb.peak(idb(j))];
        plot(or(1),or(2),'ro','buttondownfcn',{@HitPoint, j});
        hold on;
        plot(or(1)+180,or(2),'o','buttondownfcn',{@HitPoint, j});
        plot(or(1),or(2)+180,'o','buttondownfcn',{@HitPoint, j});
        plot(or(1)+180,or(2)+180,'o','buttondownfcn',{@HitPoint, j});
    end
    plot([0 360],[0 360],'-');
    axis('tight');
end


function HitPoint(a,b,cell)

fprintf('Cell%d\n',cell);

