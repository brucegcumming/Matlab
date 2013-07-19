function chspk = SetProbesToUse(DATA, ispk, distance)

nprobepc = round(distance);
if  length(ispk) > 1
    chspk = min(ispk)-nprobepc:max(ispk)+nprobepc;
elseif distance > 0 && isfield(DATA.ArrayConfig,'X')
    dp = DATA.ArrayConfig.X + i * DATA.ArrayConfig.Y;
    d = abs(dp(ispk(1)) - dp);
    chspk = find(d <= distance);
    id = find(chspk == ispk(1));
    if id ~= 2
        chspk(id) = chspk(2);
        chspk(2) = ispk(1);
    end
else
    chspk = [ispk(1)-nprobepc:ispk(1)+nprobepc];
end
chspk = chspk(chspk > 0 & chspk <= DATA.allnprobes);