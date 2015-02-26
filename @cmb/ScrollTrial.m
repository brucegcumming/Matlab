function ScrollTrial(src, evnt)

DATA = GetDataFromFig(src);

if src ~= gcf
return;
end

DATA = GetDataFromFig(src);
cmb.PlayOneTrial(DATA,src,sign(evnt.VerticalScrollCount));


