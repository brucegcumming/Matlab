function Expts = Expt2Blocks(E, varargin)
%Expts = Expt2Blocks(E, varargin) Convert a combined Expt into a cell array of individual Blocks

Bs = E.Header.BlockStart;
Be = E.Header.BlockStart(2:end);
Be(end+1) = E.Trials(end).Trial+1;
Tid = [E.Trials.Trial];
for j = 1:length(Bs)
    Expts{j} = E;
    id = find(Tid >= Bs(j) & Tid < Be(j));
    Expts{j}.Trials = E.Trials(id);
    Expts{j}.Header.BlockStart = Bs(j);
end