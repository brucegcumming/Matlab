function fign = NewFigure(tag)
fign = findobj('Tag',tag);
if ~isempty(fign)
  fign = figure('Tag',tag);
else
  figure(fign);
end

