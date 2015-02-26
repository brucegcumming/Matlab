function ReWriteProbeSpikes(dpath, varargin)

checkonly = 0;
checkfirst = 0;
matchpat = [];
j = 1;
while j<= length(varargin)
    if strncmpi(varargin{j},'check',4)
        checkonly = 1;
    elseif strncmpi(varargin{j},'match',4)
        j = j+1;
        matchpat = varargin{j};
    elseif strncmpi(varargin{j},'nocheck',4)
        checkfirst = 1;
    end
    j = j+1;
end
d = dir(dpath);
for j = 1:length(d)
    if regexp(d(j).name,'.p[0-9][0-9,t]*.mat') %% d(j).bytes < 525000000 & np < 2
        if isempty(matchpat) || length(regexp(d(j).name,matchpat))
        fprintf('loading %s\n',d(j).name)
        fname = [dpath '/' d(j).name];
        load(fname);
        a = who('-regexp','Ch*');
        if ~isempty(a) %is data in file
            err = eval([ 'size(' a{1} '.times,1) > size(' a{1} '.values,1)']);
            if checkonly && err
                fprintf('Missing values in %s\n',d(j).name);
            else
            tm=eval(['max(' a{1} '.times)']);
            if checkfirst == 0 || err
                fprintf('Saving %s to %s. Maxt %.2f\n',a{:},fname,tm);
                save(fname,a{:});
            end
            end
            clear(a{:})
        end
        end
    end
end