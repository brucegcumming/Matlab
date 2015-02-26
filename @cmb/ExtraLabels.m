function xc = ExtraLabels(Trial)
xc = '';
if isfield(Trial,'uStim') && Trial.uStim
xc = [xc 'Ustm '];
end

