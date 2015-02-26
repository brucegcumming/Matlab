function [cells, probes] = CellsSelected(DATA)
fign = findobj('Tag',DATA.tag.celllist,'Type','Figure');

np = length(DATA.probelist);
for j = 1:DATA.state.listlen
cluster(j) = 1;
it = findobj(fign,'Tag',sprintf('CellProbeAdd%d',j));
if ~isempty(it)
probes(j) = get(it(1),'value');
if probes(j) > np
probes(j) = probes(j)-np;
cluster(j) = 2;
end
end

it = findobj(fign,'Tag',sprintf('CellNumber%d',j));
if ~isempty(it)
cells(j) = get(it(1),'value')-1;
end
end
id = find(cells > 0);
cells = cells(id);
probes = probes(id);

