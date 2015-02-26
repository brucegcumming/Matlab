function KeyReleased(a, ks)
global mousept;

mousept.mode;
if strmatch(ks.Key,'delete') & mousept.mode == 5
elseif ks.Key == 'r'
cmb.ClassifySpikes(mousept, a);
end

