CheckSpikeChForValues

a = who('-regexp','Ch*');
if ~isempty(a) %is data in file
  err = eval([ 'size(' a{1} '.times,1) > size(' a{1} '.values,1)']);
  if err
    questdlg(
