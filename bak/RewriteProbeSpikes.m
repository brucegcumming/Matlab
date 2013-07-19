function ReWriteProbeSpikes(dpath)

d = dir(dpath);
for j = 1:length(d)
    if regexp(d(j).name,'.p[0-9][0-9,t]*.mat') %% d(j).bytes < 525000000 & np < 2
        fprintf('loading %s\n',d(j).name)
        fname = [dpath '/' d(j).name];
        load(fname);
        a = who('-regexp','Ch*');
        if ~isempty(a) %is data in file
            tm=eval(['max(' a{1} '.times)']);
            fprintf('Saving %s to %s. Maxt %.2f\n',a{:},fname,tm);   
            save(fname,a{:});
            clear(a{:});
        end
    end
end