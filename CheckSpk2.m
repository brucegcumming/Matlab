function details = CheckSpk2(name)

details = [];
SpkDefs;
stimch = 'Ch8';
load(name);

    vars = who;
    for j = 1:length(vars)
        if ~isempty(regexp(vars{j},'Ch[0-9][0-9]*'))
            eval(['ch = ' vars{j} ';']);
            if strncmpi(ch.title,'Spike',5)
                np = np+1;  
                probe = sscanf(ch.title,'Spike %d');
                if isempty(probe)
                    probes(np).probe = sscanf(vars{j},'Ch%d');
                    probes(np).probech = sscanf(vars{j},'Ch%d');
                else
               probes(np).probe = probe;
                end
                probes(np).var = vars{j};
            elseif strncmpi(ch.title,'StimOn',6)
                stimch = vars{j};
                stimlvl = ch;
            elseif strncmpi(ch.title,'VTR',3)
                framech = vars{j};
                fprintf('Frames in %s\n',vars{j});
            elseif strncmpi(ch.title,'Mains',5)
                mainsname = vars{j};
                mainsch = ch;
            end
        end
    end
    
out.id = find((Ch30.codes(:,1) == WURTZOK | Ch30.codes(:,1) == WURTZOKW) & Ch30.codes(:,4) == 2);
in.id = find(Ch31.codes(:,1) == WURTZOK);
out.times = Ch30.times(out.id);
in.times = Ch31.times(in.id);
rw.id = strmatch('rw',Ch30.text);
cancel.id = find(Ch31.codes(:,1) == CANCELEXPT);
rw.times = sort([Ch30.times(rw.id); Ch31.times(cancel.id)]);
bad.id = find(Ch30.codes(:,1) == BADFIX);
bad.times = Ch30.times(bad.id);
off.id = eval(['find(' stimch '.level == 0)']);
off.times = eval([stimch '.times(off.id)']);
for j = 1:length(out.id)
    id = find(rw.times >= out.times(j));
    diffs(j) = rw.times(id(1)) - out.times(j);
end
[a,b] = max(diffs);
fprintf('Max diff (%.2f) at %.2f\n',a,out.times(b));
diffs = [];
for j = 1:length(bad.id)
    id = find(off.times >= bad.times(j));
    if isempty(id)
        diffs(j) = 0;
    else
    diffs(j) = off.times(id(1)) - bad.times(j);
    end
end
[a,b] = max(diffs);
fprintf('Max (Off-Bad)(%.2f) at %.2f\n',a,bad.times(b));
hist(diffs);
id = strmatch('Man',Ch30.text);
if length(id)
    details.manstore = id;
    disp(Ch30.text(id,:));
end