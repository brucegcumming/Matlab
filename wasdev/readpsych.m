function Data = readpsych(file)


expno = 1;
stimcount = 1;
infid = fopen(file,'r');
while ~feof(infid)
  scrap = fgets(infid);
  if(scrap == -1)
    return;
  end
  resps = sscanf(scrap,'%f');
  if ~isempty(resps) %then the line begins with a digit
    Data(expno,subid).resps(stimcount,:) = resps;
    stimcount = stimcount+1;
  elseif(strncmp(scrap,'Experi',6))
    in = sscanf(scrap, 'Experiment %d %*s %*s %d')
    expno = in(1)+1;
    subid = in(2)+1;
    Data(expno,subid).label = sscanf(scrap, 'Experiment %*d %s %*s %*d');

    Data(expno,subid).sublabel = sscanf(scrap,'Experiment %*d %*s %*s %*d: %s');
    stimcount = 1;
  elseif(strncmp(scrap,'Fit',3))
    resps = sscanf(scrap,'Fit: %f %f');
    Data(expno,subid).fitsd = resps(2);
    Data(expno,subid).fitmean = resps(1);
  end
end

