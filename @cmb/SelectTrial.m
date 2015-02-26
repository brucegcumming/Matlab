function SelectTrial(a, b)
DATA = GetDataFromFig(a);
c = get(findobj(a,'Tag','ChooseTrial'),'value');
cmb.PlayOneTrial(DATA,c,0);

