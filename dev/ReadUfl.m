function res =ReadUfl(name)
%res =ReadUfl(name)
% read a .ufl file and find all bad trials.,

res.name = name;
fid = fopen(name,'r');
a = textscan(fid,'%s','delimiter','\n');
txt = a{1};
nb = 0;
ng=0;
stimid = 0;
for j = 1:length(txt)
    if strfind(txt{j},'bB')
        nb = nb+1;
        ts = sscanf(txt{j},'stim %d-%d');
        res.baddur(nb) = diff(ts);
        res.badstart(nb) = ts(1);
        res.badids(nb) = stimid;
    elseif strfind(txt{j},'bG')
        ng = ng+1;
        ts = sscanf(txt{j},'stim %d-%d');
        res.gooddur(ng) = diff(ts);
        res.goodstart(ng) = ts(1);
        res.goodids(ng) = stimid;
    elseif strncmp('props',txt{j},5) %props line appears before stim result line in ufl file
        a = strfind(txt{j},' id');
        stimid = sscanf(txt{j}(a+3:end),'%d');
    end
end