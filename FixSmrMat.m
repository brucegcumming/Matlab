function FixSmrMat(name, varargin)
%FixSmrMat(name, ..) Fix errors in a .mat file made by Spike2.
%Fixes
%If a stimon digital event seems to be missing, create a ".Stimon" file to
%correct this.


load(name);
Events = Ch30;
Text = Ch31;

vars = who('Ch*');
    for j = 1:length(vars)
        if ~isempty(regexp(vars{j},'Ch[0-9][0-9]*'))
            eval(['ch = ' vars{j} ';']);
            chn = sscanf(vars{j},'Ch%d');
            if chn > 400  %a memory/extra offline channel
            elseif strncmpi(ch.title,'SpikeO',6) && ignoreSpikeO;
            elseif strncmpi(ch.title,'Spike',5) && state.nospikes == 0
                np = np+1;  
                   if strncmpi(ch.title,'SpikeO',6)
                       probe = sscanf(ch.title,'SpikeO%d');
                       Oprobe = Oprobe+1;
                   else
                       probe = sscanf(ch.title,'Spike %d');
                   end
                if isempty(probe)
                    probes(np).probe = sscanf(vars{j},'Ch%d');
                    probes(np).probech = sscanf(vars{j},'Ch%d');
                else
                    probes(np).probe = probe;
                end
                probes(np).var = vars{j};
                probes(np).traces = ch.traces;
                probes(np).source = 1;


            elseif strncmpi(ch.title,'4Trode',5)
                np = np+1;  
                probe = sscanf(ch.title,'4Trode%d');
                if isempty(probe)
                    probes(np).probe = sscanf(vars{j},'Ch%d');
                else
               probes(np).probe = probe;
                end
                probes(np).var = vars{j};
                if probe == setprobe(1)
                    Chspk = ch;
                    spkch = 'Chspk';
                end
                probes(np).traces = ch.traces;
                probes(np).source =1;

            elseif strncmpi(ch.title,'uStimMk',7)
                ustimmarkame = vars{j};
                ustimmarkch = ch;
            elseif strncmpi(ch.title,'uStim',5)
                np = np+1;  
                probe = 100;
                UstimV = ch;
                if isempty(probe)
                    probes(np).probe = sscanf(vars{j},'Ch%d');
                    probes(np).probech = sscanf(vars{j},'Ch%d');
                else
                    probes(np).probe = probe;
                end
                probes(np).source = 0;
                probes(np).traces = ch.traces;
                probes(np).var = vars{j};
                
            elseif strncmpi(ch.title,'StimOn',6)
                stimch = vars{j};
                stimlvl = ch;
            elseif strncmpi(ch.title,'StimChan',8) % stim change detector
                dstimch = vars{j};
                stimchange = ch;
            elseif strncmpi(ch.title,'VTR',3)
                framechname = vars{j};
                framech = ch;
   %             fprintf('Frames in %s\n',vars{j});
            elseif strncmpi(ch.title,'Mains',5)
                mainsname = vars{j};
                mainsch = ch;
            elseif strncmpi(ch.title,'DigMark',7)
                ustimmarkame = vars{j};
                ustimmarkch = ch;
            end
        end
    end

    vars = who('Ch*');
    for j = 1:length(vars)
        if ~isempty(regexp(vars{j},'Ch[0-9][0-9]*'))
            eval(['ch = ' vars{j} ';']);
            chn = sscanf(vars{j},'Ch%d');
            if chn > 400  %a memory/extra offline channel
            elseif strncmpi(ch.title,'SpikeO',6) && ignoreSpikeO;
            elseif strncmpi(ch.title,'Spike',5) && state.nospikes == 0
                np = np+1;  
                   if strncmpi(ch.title,'SpikeO',6)
                       probe = sscanf(ch.title,'SpikeO%d');
                       Oprobe = Oprobe+1;
                   else
                       probe = sscanf(ch.title,'Spike %d');
                   end
                if isempty(probe)
                    probes(np).probe = sscanf(vars{j},'Ch%d');
                    probes(np).probech = sscanf(vars{j},'Ch%d');
                else
                    probes(np).probe = probe;
                end
                probes(np).var = vars{j};
                probes(np).traces = ch.traces;
                probes(np).source = 1;
            elseif strncmpi(ch.title,'4Trode',5)
                np = np+1;  
                probe = sscanf(ch.title,'4Trode%d');
                if isempty(probe)
                    probes(np).probe = sscanf(vars{j},'Ch%d');
                else
               probes(np).probe = probe;
                end
                probes(np).var = vars{j};
                if probe == setprobe(1)
                    Chspk = ch;
                    spkch = 'Chspk';
                end
                probes(np).traces = ch.traces;
                probes(np).source =1;

            elseif strncmpi(ch.title,'uStimMk',7)
                ustimmarkame = vars{j};
                ustimmarkch = ch;
            elseif strncmpi(ch.title,'uStim',5)
                np = np+1;  
                probe = 100;
                UstimV = ch;
                if isempty(probe)
                    probes(np).probe = sscanf(vars{j},'Ch%d');
                    probes(np).probech = sscanf(vars{j},'Ch%d');
                else
                    probes(np).probe = probe;
                end
                probes(np).source = 0;
                probes(np).traces = ch.traces;
                probes(np).var = vars{j};
                
            elseif strncmpi(ch.title,'StimOn',6)
                stimch = vars{j};
                stimlvl = ch;
            elseif strncmpi(ch.title,'StimChan',8) % stim change detector
                dstimch = vars{j};
                stimchange = ch;
            elseif strncmpi(ch.title,'VTR',3)
                framechname = vars{j};
                framech = ch;
   %             fprintf('Frames in %s\n',vars{j});
            elseif strncmpi(ch.title,'Mains',5)
                mainsname = vars{j};
                mainsch = ch;
            elseif strncmpi(ch.title,'DigMark',7)
                ustimmarkame = vars{j};
                ustimmarkch = ch;
            end
        end
    end

fsid = find(Events.codes(:,1) ==5 & Events.codes(:,4) ==1); %FRAMESIGNAL and Storage on
fsid = find(Events.codes(:,1) ==5); %FRAMESIGNAL and Storage on
bsid = find(stimlvl.level ==1);
esid = find(stimlvl.level ==0);
if stimlvl.level(1) == 0
    esid = esid(2:end);
end
nd =length(fsid) > length(bsid); 
fstimes = Events.times(fsid);
bstimes = stimlvl.times(bsid);
estimes = stimlvl.times(esid);
savefix = 0; %have GUI popup for this
if nd > 0
    fixfile = strrep(name,'.mat','Stimon');
    for j = 0:length(nd)
        k = nd-j;
        diffs(j+1,:) = Events.times(fsid(1+j:end-k))-stimlvl.times(bsid);
    end
    [a,b] = min(abs(diffs));
    t = find(diff(b) > 0)+1;
    dur = prctile(estimes-bstimes,90); 
    for j = 1:length(t)
        result = 0;
        if result == 0
            new = [fstimes(t) - 0.01 fstimes(t)+0.1]; %bad fix 
        else
            new = [fstimes(t) - 0.01 fstimes(t)-0.01 + dur];
        end
        s = sprintf('Missing StimOn at %.2f.  Result was %d. Adding DIO %.3f,%.3f',fstimes(t),result,new(1),new(2));
        fprintf('%s\n',s);
        newstims{j} = new;
    end
    if savefix
        BackupFile(fixfile);        
    end
    plot(diffs);
end
