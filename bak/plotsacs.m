%
%Add option to show stimulus box.circle
%add flag for X/Y, O/P and O-time plots.
function ro = plotsacs(sacs)

for j = 1:length(sacs)
  ro(j) = rotatesac(sacs(j))
  if ro(j).shift < 0
	plot(sacs(j).h,sacs(j).v,'b');
	hold on;
	plot(ro(j).spt(1),ro(j).spt(2),'ob');
  else
	plot(sacs(j).h,sacs(j).v,'r');
	hold on;
	plot(ro(j).spt(1),ro(j).spt(2),'or');
  end
end


