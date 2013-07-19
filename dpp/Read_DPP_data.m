function DPP = Read_DPP_data(trlname, forceread)

forcefile = trlname;
if exist(forcefile)
  trlname = forcefile;
end


maname = strrep(trlname,'.trl','.mma');
matname = strrep(trlname,'.trl','.mat');

if ~ exist(matname) | forceread == 1
  cmd = sprintf('trl2mat %s',trlname);
  unix(cmd);

infid = fopen(maname,'r');
title = fgets(infid);
data = [];
strings = {};
count = 1;
while ~feof(infid)
    scrap = fgets(infid);
    if(scrap ~= -1)
    trialno = sscanf(scrap,'%f');
    if ~isempty(trialno) %then the line begins with a digit
      ldata = sscanf(scrap,'%f');
      data = [data, ldata];
      count = count +1;
    else
      strings = [strings {scrap}];
    end
    end
end
fclose(infid);

DPP.title = title
DPP.data = data;
DPP.strings = strings;
if findstr(trlname,'rds')
  DPP.mainstim = 2;
elseif findstr(trlname,'rls')
  DPP.mainstim = 15;
elseif findstr(trlname,'2grating')
  DPP.mainstim = 10;
elseif findstr(trlname,'grating')
  DPP.mainstim = 3;
end

save(matname,'DPP')
cmd = sprintf('rm %s',maname);
unix(cmd);
else
    load(matname)
end
DPP.name = matname;
