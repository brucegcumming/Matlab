function cp = ExptCP(Expt)
cp = NaN;
if ~isfield(Expt,'Trials');
    return;
end
T = Expt.Trials;

if ~isfield(T,'RespDir')
    return;
end
if isfield(T,'ob')
    aid = find([T.ob] > 120 & [T.RespDir] == -1);
    bid = find([T.ob] > 120 & [T.RespDir]  == 1);
elseif isfield(T,'psyv')
    aid = find([T.psyv] == 0 & [T.RespDir] == -1);
    bid = find([T.psyv] == 0 & [T.RespDir]  == 1);
elseif isfield(T,'Dc')
    aid = find([T.Dc] == 0 & [T.RespDir] == -1);
    bid = find([T.Dc] == 0 & [T.RespDir]  == 1);
elseif isfield(T,'dx')
    aid = find([T.dx] == 0 & [T.RespDir] == -1);
    bid = find([T.dx] == 0 & [T.RespDir]  == 1);
elseif isfield(T,'pR')
    aid = find([T.pR] == 0 & [T.RespDir] == -1);
    bid = find([T.pR] == 0 & [T.RespDir]  == 1);
elseif isfield(T,'signal')
    aid = find([T.signal] == 0 & [T.RespDir] == -1);
    bid = find([T.signal] == 0 & [T.RespDir]  == 1);
elseif isfield(T,Expt.Stimvals.et)
    aid = find([T.(Expt.Stimvals.et)] == 0 & [T.RespDir] == -1);
    bid = find([T.(Expt.Stimvals.et)] == 0 & [T.RespDir]  == 1);
else
    aid = [];
    bid = [];
end
cp = CalcCP([T(aid).count],[T(bid).count]);

