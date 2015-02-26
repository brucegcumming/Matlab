function [result, details] = CheckSpkMat(varargin)

% CheckSpkMat(varargin)  checks contents of matlab variables from Spike3
% matlab files.
% CheckSpkMat(name,.....) loads file named first
% CheckSpkMat(Chx, Chy ....) Checks the matlab variables that are passed.
% CheckSpkMat(name, 'bsextra') Finds "orhhan" stimOn/Off pulses. E.g.
% lemM279.29
checkmode = 1;
result = 1;
details.errors = {};
filename = [];
fixerrs = 0;
missing = [];

if iscellstr(varargin{1})
    clear details;
    list = varargin{1};
    for j = 1:length(list)
        [result(j),details{j}] = CheckSpkMat(list{j},varargin{2:end});
        str = '';
        if isfield(details{j},'msg')
            str = details{j}.msg;
        end
        if result(j) ==1
            fprintf('%s OK%s\n',list{j},str);
        end
    end
    return;
end
j = 1;
while j <= length(varargin)
    v = varargin{j};
    if isstruct(v)
        if isfield(v,'title')
            if strfind(v.title,'Spike')
                CSpk = v;
            elseif strfind(v.title,'StimOn')
                CStim = v;
            elseif strfind(v.title,'Serial')
                CSer = v;
            elseif strfind(v.title,'Mains')
                CMains = v;
            elseif strfind(v.title,'VTR')
                CVTR = v;
            end
        end
    elseif ischar(varargin{j}) & exist(varargin{j},'file')
        filename = varargin{j};
        load(varargin{j});
        
        f = who('Ch*');
        for k = 1:length(f)
            t = eval([f{k} '.title']);
            if strfind(t,'Mains')
                CMains = eval(f{k});
            end
            if strfind(t,'StimOn')
                CStim = eval(f{k});
                StimChan = f{k};
            end
            if strfind(t,'Keyboard')
                Events = eval(f{k});
            end
            if strfind(t,'Serial')
                CSer = eval(f{k});
            end
        end
    elseif strncmpi(varargin{j},'bsextra',5)
        checkmode = 'bsextra';
    elseif strncmpi(varargin{j},'fix',3)
        fixerrs = 1;
    elseif strncmpi(varargin{j},'mains',5)
        checkmode = 2;
    elseif strncmpi(varargin{j},'vtr',3)
        checkmode = 3;
    end
    j = j+1;
end
 
if ~exist('CStim', 'var') && ischar(varargin{1})
    fprintf('Cant read %s\n',varargin{1});
    return;
end

if strcmp(checkmode,'bsextra')
    bsid = find(CStim.level == 1);
    esid = find(CStim.level == 0);
    bstimes = CStim.times(bsid);
    estimes = CStim.times(esid);
    if estimes(1) < bstimes(1)
        estimes = estimes(2:end);
        esid = esid(2:end);
    end
    dur = estimes-bstimes;
    fsid = find(Events.codes(:,1) ==5);
    for j = 1:length(bstimes)
        [a,b] = min(abs(Events.times(fsid)-bstimes(j)));
        bsdiffs(j) = Events.times(fsid(b))-bstimes(j);
    end
    gap = bstimes(2:end)-estimes(1:end-1);
    id = find(dur(2:end) < 0.02 & gap < 0.3);
    for j = length(id):-1:1
        if bsdiffs(id(j)) < 0 && (bstimes(id(j)) - bstimes(id(j)+1)) > bsdiffs(id(j));
            id(j) = [];
        end
    end
    for j = 1:length(id)
        [diffs(j),b] = min(abs(Events.times(fsid)-bstimes(id(j)+1)));
        diffs(j) = Events.times(fsid(b))-bstimes(id(j)+1);
    end
% May need to add a check that Events is really missing. 
% ie. look for Event prior to nearest and check that matches previous
% bstimes.  But foe lemM279.29, this ie enough.
    if length(id)
        ngood = length(bstimes)-length(id);
        if length(fsid) < ngood
            if bstimes(1) < Events.times(fsid(1))-1
                id = [0; id];
                ngood = ngood-1;
            end
            if length(fsid < ngood)
                xid = find(dur(2:end) < 0.02 & gap > 2 & abs(bsdiffs(2:end)') > 1);
                if length(fsid)+length(xid) == ngood
                    fprintf('Adding %d isolated On/Off. Min delay %.1f\n',length(xid),min(abs(bsdiffs(xid+1))));
                    id = [id; xid];
                elseif length(fsid)+length(xid) < ngood %still missing some
                    id = [id; xid];
                    for j = 1:length(bstimes)
                        missing(j) = sum(Events.times(fsid) < bstimes(j)+0.5) - j + sum(bstimes(id+1) <= bstimes(j));
                    end
                elseif length(fsid)+length(xid) > ngood %still missing some
                    min(abs(gaps(xid)));
                end
            end
        end
        ngood = length(bstimes)-length(id);
        x = [];
        y = [];
        for j = 1:length(bstimes)
            x = [x bstimes(j) bstimes(j) estimes(j) estimes(j)];
            y = [y 0 1 1 0];
        end
        GetFigure(['SMR.MatFile' filename]);
        hold off;
        plot(x,y);
        for j = 1:length(fsid)
            line([Events.times(fsid(j)) Events.times(fsid(j))],[0 1.1],'color','r');
        end
        x = bstimes(id+1);
        hold on;
        plot(x,ones(length(x)).*1.2,'kx');
        if ~isempty(missing)
            plot(bstimes,missing);
        end

        details.errors{end+1} = sprintf('%s Found %d Orhpan StimON Leaves %d for %d FRAME SIGNAL\n',filename,length(id),ngood,length(fsid)); 
        cprintf('error',details.errors{end})
        Fix.excludeStimOn = union(bsid(id+1),esid(id+1)); %Record fact it was fixed
        esid = setdiff(esid,esid(id+1));
        bsid = setdiff(bsid,bsid(id+1));
        goodid = union(bsid,esid);
        CStim.level = CStim.level(goodid);
        CStim.times = CStim.times(goodid);
        CStim.length = length(goodid);
        if fixerrs
            eval([StimChan '= CStim']);            
            save(filename,'Ch*','Fix');
        end
    elseif exist('Fix','var')
        details.msg = 'Previously Fixed';
    end
elseif checkmode == 1
    bstimes = CStim.times(find(CStim.level == 1));
    estimes = CStim.times(find(CStim.level == 0));
    sid = strmatch('bss',CSer.text);
    stimes= CSer.times(sid);
    eid = strmatch('ess',CSer.text);
    etimes= CSer.times(eid);
    for j = 2:length(etimes)-1;
        id = find(estimes < etimes(j)+0.015);
        delays(j,3) = etimes(j)-estimes(id(end));
        delays(j,4) = etimes(j)-estimes(id(end)+1);

        id = find(stimes < etimes(j));
        k = id(end);
        id = find(bstimes < stimes(k)+0.015);
        delays(j,1) = stimes(k)-bstimes(id(end));
        
        delays(j,2) = stimes(k)-bstimes(id(end)+1);
        starts(j) = stimes(k);
    end
    plot(delays(:,1),delays(:,2),'o')
    details.delays = delays;   
    details.starts = starts;
elseif checkmode == 2 & exist('CMains')
    [a, b] = min(diff(CMains.times));
    if a < 0.01
        fprintf('%s has short mains interval (%.6f) at %.2f\n',filename,a,CMains.times(b));
        result = 0;
        id = find(diff(CMains.times) < 0.01);
        details.times = CMains.times(id);
    end
    details.minmains = a;
elseif checkmode == 3 & exist('CVTR')
    [a, b] = min(diff(CVTR.times));
    if a < 0.005
        fprintf('%s has short VTR interval (%.6f) at %.2f\n',filename,a,CVTR.times(b));
        result = 0;
        id = find(diff(CVTR.times) < 0.01);
        details.times = CVTR.times(id);
    end
    details.minvtr = a;
else
    result = 1;
end

