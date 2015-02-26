function rsacc = plotsac(sacc)
%superceded by rotatesac
l = length(sacc.v);
rsacc.spt = [mean(sacc.h(1:10)) mean(sacc.v(1:10))]
rsacc.ept = [mean(sacc.h(l-10:l)) mean(sacc.v(l-10:l))]
rsacc.dir = atan2(rsacc.ept(2) - rsacc.spt(2),rsacc.ept(1)- ...
		  rsacc.spt(1));
dir = rsacc.dir;

rsacc.o = sacc.v -rsacc.spt(2) .* sin(dir) + sacc.h-rsacc.spt(1) .* cos(dir);
rsacc.p = sacc.h .* sin(dir) - sacc.v .* cos(dir);