function [Expts, nlines] = ReadOnlineTxt(file, DATA, varargin)
%read online data with just spike times in text format
%so that combine.m can replace binoc online plots

setprobe= 5;
setcl = 2;
SpkDefs;
nex = 0;
nt = 0;
ntall = 0;
nlines = -1;
if nargin == 1 | ~isfield(DATA,'linesread')
    nskip = 0;
else
    nskip = DATA.linesread;
    nex = length(DATA.Expts);
    if isfield(DATA.Expts{nex},'Trials')
        nt = length(DATA.Expts{nex}.Trials);
        ntall = DATA.Expts{nex}.Trials(nt).Trial;
    else
        nt = 0;
    end
    Expts = DATA.Expts;
    et = Expts{nex}.Stimvals.et;
    eb = Expts{nex}.Stimvals.e2;
    ec = Expts{nex}.Stimvals.e3;
end

fid = fopen(file);
lines = textscan(fid,'%s','delimiter','\n','HeaderLines',nskip);
lines = lines{1};
nlines = length(lines) + nskip;
expline = nlines;
emptyid = [];
for j = 1:length(lines)
    if isempty(lines{j})
        emptyid = [emptyid j];
    elseif lines{j}(1) == 'E'
        expline = j + nskip;
        if nex == 0 || isfield(Expts{nex},'Trials')
            nex = nex+1;
        end
        nt = 0;
        expts = textscan(lines{j},'%*1s%s%s%s%s%f%f');
        et = expts{1}{1};
        eb = expts{2}{1};
        ec = expts{3}{1};
        Expts{nex}.Stimvals.et = et;
        Expts{nex}.Stimvals.e2 = eb;
        Expts{nex}.Stimvals.e3 = ec;
        Expts{nex}.Stimvals.wi = expts{5}(1);
        Expts{nex}.Stimvals.hi = expts{6}(1);
        Expts{nex}.Stimvals.st = strmatch(expts{4}{1},stimnames,'exact')-1;
        Header.expname = [expts{4}{1} '.' et 'X' eb];
        Header.Name = [file];
        Header.rc = 0;
        Header.psych = 0;
        Header.trange = [0 0];
        Expts{nex}.Header = Header;
        
    elseif lines{j}(1) == 'S' && nex
        nt = nt+1;
        ntall = ntall+1;
        stim = textscan(lines{j},'%*1s%f%f%f%f%f%f%d');
        Expts{nex}.Trials(nt).Start = stim{1}.*10000;
        Expts{nex}.Trials(nt).End = stim{2}.*10000;
        if Expts{nex}.Header.trange(1) == 0
            Expts{nex}.Header.trange(1) = stim{1} .*10000;
        end
        Expts{nex}.Header.trange(2) = stim{2} .*10000; 
        Expts{nex}.Trials(nt).(et) = stim{3};
        Expts{nex}.Trials(nt).(eb) = stim{4};
        Expts{nex}.Trials(nt).(ec) = stim{5};
        Expts{nex}.Trials(nt).extra = stim{6};
        Expts{nex}.Trials(nt).id = stim{7};
        Expts{nex}.Trials(nt).OptionCode = '';
        Expts{nex}.Trials(nt).Trial = ntall;
        Expts{nex}.gui.classified = 0;
        Expts{nex}.gui.counted = 0;
        Expts{nex}.gui.clustertype = 0;
    elseif lines{j}(1) == 'R' && nex
        spikes = textscan(lines{j},'%*1s%1s%f%*[^\n]');
        if spikes{1}{1} == 'G' 
           resp = spikes{2};
        end
        if spikes{1}{1} == 'W' 
           resp = ~spikes{2};
        end
        Expts{nex}.Trials(nt).RespDir = (0.5 - resp) * 2;
    elseif lines{j}(1) == 'C' && nex
        spikes = textscan(lines{j},'%*1s%f%*[^\n]');
        probe = floor(spikes{1});
        cl = round(mod(spikes{1},1) .* 10);
        if length(lines{j}) > 5
            spktimes = textscan(lines{j}(5:end),'%f');
            Expts{nex}.Trials(nt).AllSpikes{probe,cl+1} = spktimes{1} .* 10000;
        else
            Expts{nex}.Trials(nt).AllSpikes{probe,cl+1} = [];
        end
            
     end
end

fclose(fid);
if ~isfield(Expts{nex},'Trials') %% latest expt might be empty
    Expts = {Expts{1:nex-1}};
    nex = nex-1;
    nlines = expline-1;
end
for j = 1:nex
    if isfield(Expts{j},'Trials')
        if ~isfield(Expts{j}.Trials,'me')
            [Expts{j}.Trials.me] = deal(0);
            id = find([Expts{j}.Trials.extra] == 4);
            [Expts{j}.Trials(id).me] = deal(-1);
            id = find([Expts{j}.Trials.extra] == 5);
            [Expts{j}.Trials(id).me] = deal(1);
        end
    end
end

function oldway(lines)

exl = strmatch('E',lines);
stl = strmatch('S',lines);
stims = textscan(char(lines{stl})','%*1s%f%f%f%f%f%f');
spkl = strmatch('C',lines);
spikes = textscan(char(lines{spkl})','%*1s%f%*[^\n]');
spktimes = textscan(char(lines{spkl})','%*s%*f%f');
