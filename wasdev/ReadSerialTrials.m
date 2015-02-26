function [E, txt] = ReadSerialTrials(name)
%[Trials, txt] = ReadSerialTrials(name) Read serial file fo AplaySpkfile
% convert a serial output file made by binoc in a Trials struct
% like ReadSerialFile but does not make expt strutures. 



%longfields catches any codes > 2, so that rest can be grabbed blindly

longfields = {'expname' 'puA' 'puF'  'USp' 'USd' 'USf' 'imve' 'exvals' 'imx' 'imy' 'bjv' ...
    'nbars' 'ijump' 'mixac' 'imi' 'Electrode' 'dxseq' 'ceseq' 'psyv'};
specialfields = {'EndExpt' 'EndStim' 'mt'};

ignorefields = {'Reopened' 'by binoc'  'NewConnection' 'bpos' 'sb-' 'sb+' 'Experimenter' 'Hemifield' 'StartDepth' ...
    'vebinoclean' 'Sa:' 'HemiSphere' 'CLOOP' 'Hemisphere' 'Adapter'};

charfields = {'et' 'st' 'e2' 'e3' 'Bs' 'lo'}; 

if ischar(name)
    txt = scanlines(name);
elseif iscellstr(name);
    txt = name;
end

t = 0;

Stimulus.Start = 0;

ts = now;

ignoreid = [];
for j = 1:length(ignorefields)
    id = find(strncmp(ignorefields{j},txt,length(ignorefields{j})));
    ignoreid = cat(1,ignoreid, id);
end
ignoreid = unique(ignoreid);
ignoreid(end+1) = length(txt)+1;
nextignore = 1;

specialid = [];
for j = 1:length(specialfields)
    id = find(strncmp(specialfields{j},txt,length(specialfields{j})));
    specialid = cat(1,specialid, id);
end
nextspecial = 1;

charid = [];
for j = 1:length(charfields)
    id = find(strncmp(charfields{j},txt,length(charfields{j})));
    charid = cat(1,charid, id);
end
nextchar = 1;

longid = [];
for j = 1:length(longfields)
    id = find(strncmp(longfields{j},txt,length(longfields{j})));
    longid = cat(1,longid, id);
end
longid = unique(longid);

allid = unique(cat(1,longid, ignoreid, specialid, charid));
allid(end+1) = length(txt)+1;

nextall = 1;
nextlong = 1;
nextignore = 1;
E.filename = name;
E.Trials = [];
mytoc(ts);
ts = now;
E.cellfields = [];
E.vectorfields = [];
for j = 1:length(txt);
    s = txt{j};
    if length(s) < 3
    elseif s(1) == '#'
        if strncmp(s,'#Prep',5)  %start of info describing next stimulus
            if t > 0
                f = fields(Stimulus);
                for k = 1:length(f)
                    E.Trials(t).(f{k}) = Stimulus.(f{k});
                end
            end
            t = t+1;
            [a, b,c] = sscanf(s,'#Prep%d %d %f');
        end
    elseif j == allid(nextall)
        nextall = nextall+1;
        if sum(j == longid)
            id = find(strncmp(s,longfields,3));
            f = longfields{id(1)};
            Stimulus.(f) = s(length(f)+1:end);
            if isfield(Stimulus,'EndExpt')
                f = longfields{id(1)};               
            end
        elseif sum(j == ignoreid)
        elseif sum(j == specialid)
            if strncmp(s,'EndStim',6)
                
            elseif strncmp(s,'EndExpt',6)
                nexp = nexp+1;
            end
        end
    elseif sum(strncmp(s,{'RG' 'RW' 'RB' 'RF'},2))
        if strncmp(s,'RG',2)
            Stimulus.result = 1;
        elseif strncmp(s,'RW',2)
            Stimulus.result = 1;
        elseif strncmp(s,'RB',2)
            Stimulus.result = 0;
        elseif strncmp(s,'RF',2)
            Stimulus.result = 1;
        end
    elseif strfind(s,'=')
        id = strfind(s,'=');
    elseif strfind(s,':')
        if isempty(strfind(s,'Reopened'))
        id = strfind(s,':');
        f = s(1:id(1)-1);
        E.cellfields.(f) = 1;
        a = sscanf(s(id(1)+1:end),'%f');
        Stimulus.(f) = a;
        end
    else
        f = s(1:2);
        Stimulus.(f) = sscanf(s(3:end),'%f');
        if length(Stimulus.(f)) > 1
            if sum(strcmp(f,{'sq' 'vs' 'fp'})) %Only want first value
                Stimulus.(f) = Stimulus.(f)(1);
            else
                E.vectorfields.(f) = 1;
            end
        end
        if isempty(Stimulus.(f))
            Stimulus.(f) = s(3:end);
            if ~ismember(j,charid)
                fprintf('%s\n',s);
            end
        end
    end
end
mytoc(ts);