function X = ReadSeeds(name, varargin)
%ReadSeeds  gets seed/id numbers from online text record to fix missing
%seeds in Expt files.   F
% Expts = ReadSeeds('Utah/jbe/G087/jbeG087.online',Expts,'match')
% puts correct Trial.id numbers int Expts that need fixing For MimicSaccade
% Expt
mkstims = [];
Expts = {};
stimdir = '/local/manstim';
j = 1;
while j <= length(varargin)
    if iscell(varargin{j}) && isfield(varargin{j}{1},'Trials')
        Expts = varargin{j};
    elseif strncmpi(varargin{j},'fixid',5)
        X = FixStimid(name);
        return;
    elseif strncmpi(varargin{j},'match',5)
        if ischar(name)
            X = ReadSeeds(name);
            X = ReadSeeds(X,'fixids');
            name = X;
        end
        [path, filename] = fileparts(name.name);
        R = ReadRLSFiles(path);
        for e = 1:length(Expts)
            if diff(minmax([Expts{e}.Trials.id])) == 0
                fprintf('Fixing Expt %d\n',e);
                [X, Expts{e}] = MatchTrials(Expts{e}, name, R);
            else
                ns = sum(diff([Expts{e}.Trials.id]) ==0);
                if ns > 0
                    fprintf('Expt %d %d Repeated Ids\n',e,ns);
                end
            end
        end
        X = Expts;
        return;
    elseif strncmpi(varargin{j},'mkstims',4)
        j = j+1;
        mkstims = varargin{j};
    elseif strncmpi(varargin{j},'stimdir',4)
        j = j+1;
        stimdir = varargin{j};
    end
    j = j+1;
end

fid = fopen(name,'r');
if fid < 0
    mycprintf('errors','Cant Read %s\n',name);
    return;
end
a = textscan(fid,'%s','delimiter','\n');
fclose(fid)
txt = a{1};
seid = find(strncmp('se',txt,2));
seof = find(strncmp('seof',txt,4));
seid = setdiff(seid,seof);
seof = find(strncmp('serange',txt,6));
seid = setdiff(seid,seof);
idid = find(strncmp('id',txt,2));
nfid = find(strncmp('Nf',txt,2));
stid = find(strncmp('O 5',txt,3));
endid = find(strncmp('O 3',txt,3));
said = find(strncmp('Sa',txt,2));
rbid = find(strncmp('RB ',txt,3));
badsa = find(strncmp('Sa:',txt,3));
said = setdiff(said,badsa);
respid = find(strncmp('R',txt,1));
id = find(strncmp('Ri',txt,2));
respid = setdiff(respid,id);
id = find(strncmp('Ro',txt,2));
respid = setdiff(respid,id);
id = find(strncmp('Rx',txt,2));
respid = setdiff(respid,id);
id = find(strncmp('Ru',txt,2));
respid = setdiff(respid,id);
id = find(strncmp('Ry',txt,2));
respid = setdiff(respid,id);
id = find(strncmp('Rw',txt,2));
respid = setdiff(respid,id);
id = find(strncmp('Rh',txt,2));
respid = setdiff(respid,id);
id = find(strncmp('Re',txt,2));
respid = setdiff(respid,id);

badsa = [];
for j = 1:length(said)
    if isempty(strfind(txt{said(j)},','))
        badsa = [badsa said(j)];
    end
end
said = setdiff(said,badsa);

stimid = 0;

for j = 1:length(stid)
    X.start(j) = sscanf(txt{stid(j)},'O 5 %f'); 
    es = find(endid > stid(j));
    if ~isempty(es)
        es = endid(es(1));
    else
        es = NaN;
    end
    eb = find(rbid > stid(j));
    if ~isempty(eb) && rbid(eb(1)) < es
        es = rbid(eb(1));
    end
    if j < length(stid)
        next = stid(j+1);
        id = seid(find(seid < next & seid > es));
        if length(id) > 1 && id(1)+1 ==  id(2)
            X.firstseeds(j) = sscanf(txt{id(1)},'se%d');
            X.lastseeds(j) = sscanf(txt{id(2)},'se%d');
        else
            id = seid(find(seid < stid(j)));
            X.firstseeds(j) = sscanf(txt{id(end)},'se%d');
            id = find(respid < next);
            s = txt{respid(id(end))};
            si = findstr(s,'se=');
            X.lastseeds(j) = sscanf(s(si:end),'se=%d');
        end
        
        if es < next
            next = es;
        end
        id = said(find(said < next));        
        if ~isempty(id)
            x = regexprep(txt{id(end)},'.*,','');
            stimid = sscanf(x,'%d');
        end
    end
    X.stimids(j) = stimid;
    id = idid(find(idid < stid(j)));
    X.ids(j) = sscanf(txt{id(end)},'id%d');;
    id = seid(find(seid < stid(j)));
    X.seeds(j) = sscanf(txt{id(end)},'se%d');;
    id = nfid(find(nfid > stid(j)));
    if isempty(id)
    X.nf(j) = NaN;
    else
        X.nf(j) = sscanf(txt{id(1)},'Nf%d');
    end
end
X.name = name;

if ~isempty(Expts)
    for e = 1:length(Expts)
        for t = 1:length(Expts{e}.Trials)
            id = find(X.ids == Expts{e}.Trials(t).id);
            fid = fopen([stimdir '/stim' num2str(t-1)],'w');
            fprintf(fid,'id%d\nse%d\nnf%d\n',X.ids(id),X.seeds(id),X.nf(id));
            fclose(fid);
        end
        fid = fopen([stimdir '/stimorder'],'w');
        for t = 1:length(Expts{e}.Trials)
            fprintf(fid,'%d ',t-1);
        end
        fprintf(fid,'\n');
        fclose(fid);
    end
end
            
if ~isempty(mkstims)
    for j = 1:length(mkstims)
        fid = fopen([stimdir '/stim' num2str(j)],'w');
        fprintf(fid,'id%d\nse%d\n',X.ids(j),X.seeds(j));
        fclose(fid);
    end
end

function X = FixStimid(X)

X.fixids = X.stimids;
id = find(X.ids > X.stimids);
X.fixids(id) = X.ids(id);
id = find(diff(X.stimids) ==0);
xid = diff(X.ids);
for j = 1:length(id)
    if xid(id(j)) == 0
        X.fixids(id(j)+1) = X.fixids(id(j))+1;
    else
        x = X.ids(id(j));
    end        
end
id = find(diff(X.fixids) ==0);
for j = 1:length(id)
end


function [M, Expt] = MatchTrials(Expt, X, R)

allR = cat(1,R{:});
    
for j = 1:length(Expt.Trials);
    t = Expt.Trials(j).Start(1)./10000;
    [a,b] = min(abs(t-X.start));
    sid = find(X.firstseeds == Expt.Trials(j).se(1));
    g = find(ismember(sid,b));
    if isempty(g)
        M.match(j) = 0;
        Expt.Trials(j).matchid = 0;
    else
        M.match(j) = sid(g);
        Expt.Trials(j).id = X.fixids(sid(g));
        Expt.Trials(j).matchid = find(allR(:,1) == X.fixids(sid(g)) & allR(:,2) == X.lastseeds(sid(g)));
    end
end

function X = ReadRLSFiles(path)

X = {};
d = mydir([path '/stims/rls*']);
for j = 1:length(d)
    txt = scanlines(d(j).name);
    idid = find(strncmp('id',txt,2));
    for t = 1:length(idid)
        X{j}(t,:) = sscanf(txt{idid(t)},'id%dse%d');
    end
end

    