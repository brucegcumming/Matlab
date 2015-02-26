function X = ReadSeeds(name, varargin)
%ReadSeeds  gets seed/id numbers from online text record to fix missing
%seeds in Expt files.   F
% Expts = ReadSeeds('Utah/jbe/G087/jbeG087.online',Expts,'match')
% puts correct Trial.id numbers int Expts that need fixing For MimicSaccade
% Expt
mkstims = [];
Expts = {};
savedx = 1;
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
id =  find(strncmp('Rh',txt,2));
respid = setdiff(respid,id);
id = find(strncmp('Re',txt,2));
respid = setdiff(respid,id);
dxid = find(strncmp('dx:',txt,3));
ceid = find(strncmp('ce:',txt,3));
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
    es = find(endid < stid(j));
    if ~isempty(es)
        es = endid(es(end)); %preceding End
    else
        es = NaN;
    end
    eb = find(rbid < stid(j));
    if ~isempty(eb) && rbid(eb(end)) > es
        es = rbid(eb(end));
    end
    id = dxid(find(dxid > es & dxid < stid(j)));
    if ~isempty(id)
        x = regexprep(txt{id(end)},'.*,','');
        dxvals = sscanf(x(4:end),'%f');
        X.dxvals{j} = dxvals;
    end
    id = ceid(find(ceid > es & ceid < stid(j)));
    if ~isempty(id)
        x = regexprep(txt{id(end)},'.*,','');
        cevals = sscanf(x(4:end),'%f');
        X.cevals{j} = cevals;
    end
    id = idid(find(idid > es & idid < stid(j)));
    if ~isempty(id)
        stimid = sscanf(txt{id(end)},'id%f');
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

if savedx
    outname = strrep(name,'online','SetStim');
    fid = fopen(outname,'w');
    for j = 1:length(X.dxvals)
        if ~isempty(X.dxvals{j})
            n = min([length(X.dxvals{j}) X.nf(j)]);
            fprintf(fid,'id%d dx:%s\n',X.ids(j),sprintf('%.3f ',X.dxvals{j}(1:n)));
            fprintf(fid,'id%d ce:%s\n',X.ids(j),sprintf('%.1f ',X.cevals{j}(1:n)));
        end
    end
    fclose(fid);
end

if ~isempty(Expts)
    
    for e = 1:length(Expts)
        for t = 1:length(Expts{e}.Trials)
        end
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

    