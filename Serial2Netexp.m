function [E, txt] = Serial2Netexp(name)
%[Trials, txt] = Serial2Netexp(name) Read serial file, and create
%file that looks like network stim file that binoc would have made...
%Saves result to name.bno. Rename this to name.bnc an rerun AplaySpkFile
%to get the changes into the .mat files


%longfields catches any codes > 2, so that rest can be grabbed blindly

longfields = {'expname' 'puA' 'puF'  'USp' 'USd' 'USf' 'imve' 'exvals' 'imx' 'imy' 'bjv' ...
    'nbars' 'ijump' 'mixac' 'imi' 'Electrode' 'dxseq' 'ceseq' 'psyv'};
specialfields = {'EndExpt' 'EndStim' 'mt'};

ignorefields = {'Reopened' 'by binoc'  'NewConnection' 'bpos' 'sb-' 'sb+' 'Experimenter' 'Hemifield' 'StartDepth' ...
    'vebinoclean' 'Sa:' 'HemiSphere' 'CLOOP' 'Hemisphere' 'Adapter' 'Codes:' 'Outcodes:'};

charfields = {'et' 'st' 'e2' 'e3' 'Bs' 'lo'}; 

oid = -1;
if ischar(name)
    txt = scanlines(name);
    outname = regexprep(name,'.online','.bno');
    if strcmp(outname,name);
        fprintf('Cant use same name for output\n');
        return;
    end
    oid = fopen(outname,'w');
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
E.Trials = [];
mytoc(ts);
ts = now;
E.cellfields = [];
E.vectorfields = [];
for j = 1:length(txt);
    s = txt{j};
    if length(s) < 3
    elseif s(1) == '#'
        if strncmp(s,'#du',3) || strncmp(s,'#Prep',5) 
            fprintf(oid,'%s\n',s);
        end
    elseif sum(strncmp(s,{'id' 'se' 'Nf'},2))
        fprintf(oid,'%s\n',s);
    elseif sum(strncmp(s,{'RG' 'RW' 'RB' 'RF'},2))
        fprintf(oid,'%s\n',s);        
    elseif sum(strncmp(s,{'mtFl=' 'mtFn=' 'mtFi='},4))
        fprintf(oid,'%s\n',s);        
    elseif strfind(s,':')
        fprintf(oid,'%s\n',s);        
    elseif regexp(s,'exvals[0-9]')
        fprintf(oid,'#Prep\n');
    end
end
fclose(oid);
mytoc(ts);