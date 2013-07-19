function mnames = GetTFNames(file,idstr)

mnames = {};
fcnt = 1;

suffixes = {'TF', 'TFM', 'CTF'};
grfile = strrep(file,'rds','grating');
prefixes = {'grfile', 'file'};

for k = 1:2
  for j = 1:length(suffixes)
  tffile = strrep(grfile,idstr,suffixes{j});
  if(exist(tffile,'file'))
    mnames{fcnt} = tffile;
    fcnt = fcnt + 1;
  end
end
end


