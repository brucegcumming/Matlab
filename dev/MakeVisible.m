function MakeVisible(F)


x = get(F,'Position');
x(1) = 1;
x(2) = x(4);
set(F,'Position',x);