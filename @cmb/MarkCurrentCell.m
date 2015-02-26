function DATA = MarkCurrentCell(DATA)

square = [-.5 0.5 0.5 -0.5 -0.5; -0.5 -0.5 0.5 0.5 -0.5];
hold on;
if isfield(DATA, 'markh') && ishandle(DATA.markh)
delete(DATA.markh);
end
DATA.markh = plot(square(1,:)+DATA.currentpoint(2),square(2,:)+DATA.currentpoint(1),'w-','linewidth',2);

