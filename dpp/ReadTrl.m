infid = fopen('testinput','r');
title = fgets(infid);
data = [];
strings = [];
while ~feof(infid)
    scrap = fgets(infid);
    trialno = sscanf(scrap,'%f');
    if ~isempty(trialno) %then the line begins with a digit
      ldata = sscanf(scrap,'%f');
      data = [data, ldata];
    else
      strings = [strings, scrap];
    end
end
fclose(infid);



dumean = mean(data(4,:));
dpvals = data(16,:);
dqvals = data(17,:);