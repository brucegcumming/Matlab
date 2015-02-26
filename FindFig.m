function F = FindFig(tag)

F = findobj(get(0,'children'),'flat','type','figure','tag',tag);