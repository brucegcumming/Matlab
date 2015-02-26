function myjcontrol(varargin);


F = GetFigure('Jcontrol');
delete(allchild(F));
s  = 'Start';
for j = 1:100
    str{j} = sprintf('Line%d',j);
    s = sprintf('%s\n%s',s,str{j});
end
scrollpane=jcontrol(F, 'javax.swing.JScrollPane', 'Position', [0 0 0.98 0.98]);
lst =javax.swing.JEditorPane('Text', s);
scrollpane.setViewportView(lst);

