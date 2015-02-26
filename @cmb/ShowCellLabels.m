function ShowCellLabels(DATA)
if ~isfield(DATA,'CellList')
return;
end
probe = GetProbe(DATA, DATA.currentexpt(1), DATA.probe);
if ndims(DATA.CellList) == 3 && DATA.currentexpt(1) <= size(DATA.CellList,1)
cells = squeeze(DATA.CellList(DATA.currentexpt(1), probe,:));
yl = get(gca,'ylim');
xl = get(gca,'xlim');
dh = 0;
for j = 1:length(cells)
if cells(j) > 0
h =  text(xl(2),yl(1)+dh,sprintf('Cell%d',cells(j)),...
'VerticalAlignment','Bottom','HorizontalAlignment','Right','color',DATA.spkcolor{j+1});
sz = get(h,'Extent');
dh = dh + sz(4);
end
end
end

