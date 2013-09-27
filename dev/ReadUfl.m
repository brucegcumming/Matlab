function res =ReadUfl(name)
%res =ReadUfl(name)
% read a .ufl file and find all bad trials.,

res.name = name;

nb = 0;
ng=0;
stimid = 0;

if isdir(name)
    d = mydir([name '/*.ufl']);
    nrf = 0;
    rfs = [];
    dates = [];
    for j = 1:length(d)
        [rf, date] = ReadRF(d(j).name);
        if ~isempty(regexp(d(j).name,'[M][0-9][0-9]'))
            types(j) = 1;
        else
            types(j) = 2;
        end
        if ~isempty(rf)
            rfs = [rfs rf];
            if ~isempty(date)
            dates(size(rfs,1)) = date;
            end
        end
    end
    res.rfs = rfs;
    res.dates = dates;
else
    [rf, date, txt] = ReadRF(name);
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
end

function [rfs, dates, txt] = ReadRF(name)

rfs = [];
dates = [];
nrf = 0;
txt = scanlines(name);
rid = find(strncmp('cm=rf',txt,5));
for k = 1:length(rid)
    nrf = nrf+1;
    rf{nrf} = sscanf(txt{rid(k)},'cm=rf%f,%f:%fx%f,%fdeg pe%f %f,%f fx=%f,fy=%f');
    rfs(nrf,1:length(rf{nrf})) = rf{nrf};
end
did = find(strncmp('Created',txt,7));
if ~isempty(did) && length(txt{did(end)}) > 13
    dates = datenum(txt{did(end)}(9:end));
end
