function [E, txt] = ReadSerialTrials(name)
%[Trials, txt] = ReadSerialTrials(name) Read serial file fo AplaySpkfile
% convert a serial output file made by binoc in a Trials struct
% like ReadSerialFile but does not make expt strutures. 

maxchar = 100;

%longfields catches any codes > 2, so that rest can be grabbed blindly

longfields = {'expname' 'puA' 'puF'  'USp' 'USd' 'USf' 'imve' 'exvals' 'bjv' ...
    'nbars' 'ijump' 'mixac' 'imi' 'Electrode' 'dxseq' 'ceseq' 'psyv' 'seof' 'serange'};
specialfields = {'EndExpt' 'EndStim' 'mt'};

ignorefields = {'Reopened' 'by binoc'  'NewConnection' 'bpos' 'sb-' 'sb+' 'Experimenter' 'Hemifield' 'StartDepth' ...
    'vebinoclean' 'Sa:' 'HemiSphere' 'CLOOP' 'Hemisphere' 'Adapter' 'uf'};

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
lastid = 0;
Stimulus.idoffset = 0;
lastprep = 0;
instim = 0;
nchar = 0;

% for each trial outcome line ('RG' 'RB' 'RF' 'RW')
%keep reading until the next #prep line (Start of next stim)
%Can sometimes get #prep twice with no RG in between. When this 
%happens it all belongs to the next R[GWB] line
for j = 1:length(txt);
    s = txt{j};
    if length(s) < 3
    elseif s(1) == '#'
        if sum(strncmp(s,{'#Prep' '#Resetting'},5))  && instim %start of info describing next stimulus or end of expt block
            instim = 0;
            if t > 0 && (j - lastprep) > 6 %if made from online, this line might be repeated
                if Stimulus.id < lastid
                    Stimulus.idoffset = lastid + Stimulus.idoffset;
                elseif Stimulus.id == lastid %repeat should not happen
                    fprintf('Repeated id %d on line %d. Last prep line %d \n',lastid,j,lastprep);
                end
                lastid = Stimulus.id;
                f = fields(Stimulus);
                for k = 1:length(f)
                    E.Trials(t).(f{k}) = Stimulus.(f{k});
                    if length(Stimulus.(f{k})) > 5 && ~ischar(Stimulus.(f{k}))
                        Stimulus.(f{k})  = []; %don't allow values to carry over
                    elseif sum(strncmp(f{k},{'mtFn' 'mtFl' 'mtFi'},4))
                        Stimulus.(f{k})  = []; %don't allow values to carry over                            
                    end
                end
                if isfield(Stimulus,'ce') && length(Stimulus.ce) > 20
                    Stimulus.ce = 1; %reset to ensure that ce seq doesn't persist
                end
                t = t+1;
            end
            if t == 0
                t = 1;
            end
            if isfield(Stimulus,'mtFn') && ~isempty(Stimulus.mtFn)
%                diff(Stimulus.Fi);
            end
            lastprep = j;
            [a, b, c] = sscanf(s,'#Prep%d %d %f');
        elseif strncmp(s,'#du',3)  %start of info describing next stimulus
            a = sscanf(s,'#du%f(%f:%f');
            Stimulus.duration = a(1);
            if length(a) > 1
                Stimulus.nframes = a(2);
            end
            if length(a) > 2
                Stimulus.nomdur = a(3);
            end
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
            elseif strncmp(s,'mtFl=',5)
                Stimulus.mtFl = sscanf(s(6:end),'%f');
                E.cellfields.mtFl = 1;
            elseif strncmp(s,'mtFn=',5) && length(s) > 5
                Stimulus.mtFn = sscanf(s(6:end),'%f');
                E.cellfields.mtFn = 1;
            elseif strncmp(s,'mtFi=',5)
                Stimulus.mtFi = sscanf(s(6:end),'%f');
                E.cellfields.mtFi = 1;
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
        instim = 1;
    elseif strncmp(s,'O ',2)    
    elseif strfind(s,'=')
        id = strfind(s,'=');
    elseif strfind(s,'Log:')
    elseif strncmp(s,'Remaining',8)
    elseif strfind(s,':')
        if ~isempty(strfind(s,'Reopened'))
        elseif ~isempty(strfind(s,'Resetting'))
        elseif ~isempty(strfind(s,'Run ended'))
        elseif strncmp(s,'Codes:',5)    
        elseif strncmp(s,'Outcodes:',5)    
        elseif strncmp(s,'epos:',5)    
        elseif strncmp(s,'Unrecognized',9)    
        else
        id = strfind(s,':');
        f = s(1:id(1)-1);
        if isvarname(f)  %If not a valid field name his is not a real sequence parameter
            E.cellfields.(f) = 1;
            a = sscanf(s(id(1)+1:end),'%f');
            Stimulus.(f) = a;
        else
            fprintf('Ignoring %s\n',s);
        end
        end
    else
        if sum(strncmp(s,{'imx' 'imy'},3))
            f = s(1:3);
            nc = 4;
        else
            f = s(1:2);
            nc = 3;
        end
        if strcmp(f,'ce') && isfield(Stimulus,'ce') && length(Stimulus.ce) > 20
            Stimulus.ceval = sscanf(s(3:end),'%f');
        else
            try %in case of invalid fields
            Stimulus.(f) = sscanf(s(nc:end),'%f');
            end
        end
        if isfield(Stimulus,f)
        if length(Stimulus.(f)) > 1
            if sum(strcmp(f,{'sq' 'vs' 'fp'})) %Only want first value
                Stimulus.(f) = Stimulus.(f)(1);
            else
                E.vectorfields.(f) = 1;
            end
        end
        if sum(strncmp(s,{'mtFi' 'mtFn'},4))
%                E.vectorfields.(f) = 1; %cell field in Trials.xx for idx
        end
            
        if isempty(Stimulus.(f))
            Stimulus.(f) = s(nc:end);
            if ~ismember(j,charid)  && nchar < maxchar
                fprintf('%s\n',s);
                nchar = nchar+1;
            end
        end
        end
    end
end
E.readdate = now;
mytoc(ts);