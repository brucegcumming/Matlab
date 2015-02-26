function ClearOSpike(dname, varargin)

if isdir(dname)
    d = dir(dname);
end

for j = 1:length(d)
    if regexp(d(j).name,'\.[0-9,a]*\.mat')
        clear Ch*;
        load([dname '/' d(j).name]);
        chs = who('Ch*');
        chnames = {};
        for k = 1:length(chs)
            isspk = eval(['isfield(' chs{k} ',''values'')']);
            if isspk == 0
                chnames = {chnames{:} chs{k}};
            end
        end
        if length(chnames) > 0 && length(chs) > length(chnames)
            fprintf('Re-writing %s\n',d(j).name);
        save([dname '/' d(j).name],chnames{:});
        end
    end
end