%CheckSpikeChForValues

a = who('-regexp','Ch.*');
if ~isempty(a) & isfield(a{1},'times') & isfield(a{1},'values') %is data in file
  err = eval([ 'size(' a{1} '.times,1) > size(' a{1} '.values,1)']);
  if err
    if exist('filename','var')
    questdlg(sprintf('No Values in %s',filename),a{1},'OK','OK')
    else
      questdlg(sprintf('No Values in file'),a{1},'OK','OK')
    end
  end
end
clear a;

    
