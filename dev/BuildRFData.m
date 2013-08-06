function rfs = BuildRFData(dirname, varargin)

recurse = 0;
rfs = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'find',4)
        recurse = 1;
    end
    j = j+1;
end

if recurse
    nc = 0;
    d = mydir(dirname);
    for j= 1:length(d)
        if d(j).isdir
            rf = GetRFFromDir(d(j).name);
            if ~isempty(rf)
                nc = nc+1;
                rfs{nc} = rf;
            end
        end
    end
    for j = 1:length(rfs)
    end
end


function therf = GetRFFromDir(dirname)
d = mydir([dirname '/*.ufl']);
therf = [];

if isempty(d)
    return;
end
nrf = 0;
for j = 1:length(d);
    txt = scanlines(d(j).name);
    rid = find(strncmp('cm=rf',txt,5));
    for k = 1:length(rid)
        nrf = nrf+1;
        rf{nrf} = sscanf(txt{rid(k)},'cm=rf%f,%f:%fx%f,%fdeg pe%f %f,%f fx=%f,fy=%f');
        rfs(nrf,1:length(rf{nrf})) = rf{nrf};
    end
    if ~isempty(regexp(d(j).name,'[M][0-9][0-9]'))
        types(j) = 1;
    else
        types(j) = 2;
    end
end
type = prctile(types,50);
if type == 1
    therf.electrode= 'uProbe';
else
    therf.electrode = 'Normal';
end
therf.rf = prctile(rfs,50);
therf.name = dirname;
