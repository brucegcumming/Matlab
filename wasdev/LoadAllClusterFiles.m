function LoadAllClusterDetails(name, varargin)

rewrite = 0;
readnew = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'rewrite',5)
        rewrite = 1;
    elseif strncmpi(varargin{j},'readnew',5)
        readnew = 1;
    end
    j = j+1;
end

if readnew
    d = mydir([name '/Tests/*ClusterTimesDetails.mat']);
else
    d = mydir([name '/*ClusterTimesDetails.mat']);
end

ts = now;
for j = 1:length(d)
    if strfind(d(j).name,'Details')
        xname = d(j).name;
    else
        xname = strrep(d(j).name,'ClusterTimes','ClusterTimesDetails');
    end
    if exist(xname,'file');
        load(xname);
        if rewrite
            [a,b,c] = fileparts(xname);
            outname = [name '/Tests/' b c];
            for c = 1:length(ClusterDetails)
                if isfield(ClusterDetails{c},'xy')
                    X.maxint = 32000;
                    scale = max(ClusterDetails{c}.xy(:))./X.maxint;
                    X.scale(1) = scale;
                    X.xy = int16(round(ClusterDetails{c}.xy .*scale));                    
                    scale = max(ClusterDetails{c}.t)./X.maxint;
                    X.scale(2) = scale;
                    X.t = int16(ClusterDetails{c}.t);
                    X.clst = int16(ClusterDetails{c}.clst);
                    ClusterDetails{c} = X;
                end
            end
            save(outname,'ClusterDetails');
        end
        if readnew
            for c = 1:length(ClusterDetails)
                C = ClusterDetails{c};
                if isfield(C,'xy')
                ClusterDetails{c}.xy = double(C.xy) .* C.scale(1);
                ClusterDetails{c}.t = double(C.t) .* C.scale(2);
                end
            end
        end
    end
end
mytoc(ts);
