function TestEvents(name)
tic;
load(name);
Events = Ch31;
Text = Ch30;
SpkDefs;
nt = 0;
nx = 0;
nonstore = 0;
findtrial = 0;
nomdur = 5000;
Expts = [];
starttrial = 0;
Expt.setprobe = 0; 
frametimes = [];
bstimes = [];
estimes = [];
Events.times = Events.times .* 10000;
Text.times = Text.times .* 10000;
lastend = 0;

toc;
tic;
frametimes = Ch17.times * 10000;
stimch = 'Ch36';
stimlvl = Ch36;
if exist(stimch,'var') & isfield(stimlvl,'level') & length(stimlvl.times) > 1
%    id = find(Ch8.level == 1);
% need to add an extra event at the end of each list to avoid issues with
% empty finds on the last trial. But make it long after to avoid any
% confusion with the real one;
    if isfield(stimlvl,'inverted') && stimlvl.inverted
        bstimes = stimlvl.times(stimlvl.level == 0) * 10000;
        bstimes(end+1) = bstimes(end)+10000;
        estimes = stimlvl.times(stimlvl.level == 1) * 10000;
        estimes(end+1) = estimes(end)+10000;
    else
        bstimes = stimlvl.times(stimlvl.level == 1) * 10000;
        bstimes(end+1) = bstimes(end)+10000;
        estimes = stimlvl.times(stimlvl.level == 0) * 10000;
        estimes(end+1) = estimes(end)+10000;
    end
end


opid = strmatch('op',Text.text);
for j = 1:length(opid)
    str= sscanf(Text.text(opid(j),:),'op%d');
    if ~isempty(str)
        storing(j) = str;
    else
        storing(j) = 0;
    end
end
storing = bitand(storing,STOREBIT);
%storeonoff is a list of times where storing is toggled
storeonoff = 1+find(abs(diff(storing)) > 0); 
if isempty(storeonoff)
    storeonoff = [1 1];
end
instim = 0;
inexpt = 0;
nextonoff = 1;
    storestate = storing(1);
    tonoff = Text.times(opid(storeonoff(nextonoff)));
Events.store = zeros(size(Events.times));
if ~storing(storeonoff(1)) & storing(1) %% first event is an off
    onid = find(Events.times < Text.times(opid(storeonoff(1))));
    Events.store(onid) = 1;
end

for j = 1:length(storeonoff)
    if storing(storeonoff(j))
        if length(storeonoff) > j
        onid = find(Events.times >= Text.times(opid(storeonoff(j))) ...
            & Events.times < Text.times(opid(storeonoff(j+1))));
        else
        onid = find(Events.times >= Text.times(opid(storeonoff(j))));
        end
        Events.store(onid) = 1;
    end
end
fprintf('Store Index takes %.2f\n',toc);

tic;
for j = 1:size(Events.codes,1)

% if storage turned of mid-stim, don't want to miss ENSTIM marker
    if Events.codes(j,1) ~= ENDSTIM
        storestate = Events.store(j);
    end        
    if Events.codes(j,1) == FRAMESIGNAL
        nowt = Events.times(j);
        if storestate
            if nt & Trials.Result(nt) < 0
                nt = nt;
            end
            nt = nt+1;
        Trials.Start(nt) = Events.times(j);
        Trials.End(nt) = Trials.Start(nt)+nomdur; %% just in case an online file is missing end
        if findtrial  & Trials.Start(nt) > findtrial
            findtrial = 0;
        end
        Trials.Startev(nt) = Events.times(j);
        Trials.stored(nt) = storestate;
        Trials.Result(nt) = 1;
        if nt > 1 & length(Trials.End) < nt-1
            nerr = nerr+1;
            errs{nerr} = sprintf('Missing end Trial %d (EX %.0f, start %.0f)\n',nt-1,inexpt,Trials.Start(nt-1));
            fprintf('Missing end Trial %d (EX %.0f, start %.0f)\n',nt-1,inexpt,Trials.Start(nt-1));
%            Trials.End(nt-1) = NaN;
        end
        Trials.Trial(nt) = nt+starttrial;
        instim = 1;
%the event time can be just before the stim signal since it is not delayed
%to the sync (? correct explanation. Definitely what happens
% is the times of teh StimulusON digital event markers
%frametimes are the Vsync digital markers (every frame)
% the digital step can be >20ms after the serial signal is received, e.g.
% in ruf1989 at 209.43 sec
% in lem017 at 120.6 it is nearl 80ms late. at 289.6 its 300ms late. Looks
% like we could check for immediately preceding STARTSTIM (6) being before 
% StimON to check for this

        id = find(bstimes < Events.times(j)+300);
        
        if Events.times(j) > 238170000
            bstimes(id(end));
        end
        if id
            Trials.serdelay(nt) = Events.times(j) - bstimes(id(end));
            Trials.bstimes(nt) = bstimes(id(end));
        end
        if ~isempty(id) & bstimes(id(end)) > lastend 
            if id(end) < length(bstimes) & bstimes(id(end)+1) - Trials.Start(nt) < 10 %< 1ms to next = probably early
                Trials.Start(nt) = bstimes(id(end)+1);
            elseif bstimes(id(end)) > lastend
                Trials.Start(nt) = bstimes(id(end));
            end
            if Trials.serdelay(nt) > 10000
                nt = nt;
            end
            Trials.stored(nt) = storestate;
%here, Trials.Start is the time of the Digital event marker
%If vertical retrace was recorded, set start time to next one of these
            if ~isempty(frametimes) % have VTR channel
                id = find(frametimes > Trials.Start(nt) & frametimes < Trials.Start(nt)+500);
                if ~isempty(id)
                    Trials.delay(nt) = frametimes(id(1)) - Trials.Start(nt);
                    Trials.Start(nt) = frametimes(id(1));
                    Trials.FalseStart(nt) = 0;
                else
                    Trials.FalseStart(nt) = 1;
                    Trials.delay(nt) = NaN;
                end
            else
%                Trials.Start(nt) = Events.times(j);
                Trials.FalseStart(nt) = 2;
                Trials.delay(nt) = NaN;
            end
        else
            Trials.Start(nt) = Events.times(j);
            Trials.delay(nt) = NaN;
            if isempty(id)
                Trials.FalseStart(nt) = 1;
            elseif id(end) < length(bstimes) & bstimes(id(end)+1) - Events.times(j) < 800 ...
%                    & Events.codes(j-1) == STARTSTIM ... %< 1ms to next = probably early
            nerr = nerr+1;
                errs{nerr} = sprintf('StimON at %.2f is %.1f ms late but STARTSTIM at %.2f',...
                bstimes(id(end)+1),(bstimes(id(end)+1)-Events.times(j))./10,Events.times(j-1));
            fprintf('%s\n',errs{nerr}); 
                Trials.FalseStart(nt) = 0;
                
            else
%Serial input can be very late if Spike2 got busy. Use the DIO stimon -
%this is the true start
                Trials.FalseStart(nt) = Events.times(j) - bstimes(id(end));
                Trials.Start(nt) = bstimes(id(end));
                 nerr = nerr+1; errs{nerr} = sprintf('Missing StimON at %.2f %.2f), but STARTSTIM at %.2f',Events.times(j),bstimes(id(end)),Events.times(j-1));
                fprintf('%s\n',errs{nerr});
                id = find(frametimes > Trials.Start(nt));
                if ~isempty(id) 
                    Trials.delay(nt) = frametimes(id(1)) - Trials.Start(nt);
                    Trials.Start(nt) = frametimes(id(1));
                end
            end
            Trials.stored(nt) = storestate;
        end 
        else
            nonstore = nonstore+1;
        end %if storestate
  
    elseif Events.codes(j,1) == ENDSTIM & storestate & nt
    if abs(Events.times(j) - 62740381) < 200
            Trials.End(nt) = Events.times(j);
            Trials.endelay(nt) = NaN;
    end
        if estimes
            id = find(estimes < Events.times(j)+500);
            if(id)
                Trials.TrueEnd(nt) = estimes(id(end));
                Trials.End(nt) = estimes(id(end));
                Trials.endelay(nt) = NaN;
                Trials.estimes(nt) = estimes(id(end));
  %if this is out by 400ms, probably failed to find correct end mark
                if Trials.TrueEnd(nt) < Events.times(j) - 4000  
                    fprintf('End event %.3f but marker %.3f\n',...
                        Events.times(j)./10000,Trials.TrueEnd(nt)./10000);
                    if Trials.End(nt) < Trials.Start(nt) & length(estimes) > id(end) & ...
                        estimes(id(end)+1) - Events.times(j) < 10000
                        Trials.TrueEnd(nt) = estimes(id(end));
                        Trials.End(nt) = estimes(id(end));
                        Trials.endelay(nt) = NaN;
                        Trials.estimes(nt) = estimes(id(end));
                    end
                end

                if ~isempty(frametimes) % have VTR channel
                    id = find(frametimes > Trials.End(nt));
                    if ~isempty(id) & frametimes(id(1))-Trials.TrueEnd(nt) < 500
                        Trials.endelay(nt) = Trials.End(nt) - frametimes(id(1));
                        Trials.End(nt) = frametimes(id(1));
                        Trials.TrueEnd(nt) = frametimes(id(1));
                    else
                    end
                end
            else
                Trials.TrueEnd(nt) = 0;
            end
            
        end
        Trials.End(nt) = Events.times(j);
        Trials.Result(nt) = 1;
        if (Trials.End(nt) - Trials.Start(nt)) < 1000
            instim = 0;
        end
        instim = 0;
        if Trials.TrueEnd(nt)
            lastend = Trials.TrueEnd(nt);
        else
            lastend = Trials.End(nt);
        end
    elseif Events.codes(j,1) == ENDTRIAL & storestate
        if instim
%  can't figure this out here because the BADFIX is only recorded in text,
%  not SampleKey (becuase this is send from Spike2, not received by, and
%  setting codes for sample keys is such a pain. But maybe should make all
%  of these events with code2 set to indicate it is from Spike2?
% Seems like this happens when fixation is broken just BEFORE stimulus on,
% but Spike2 has not registered this yet e.g. ruf2000 at 8867.9
%            fprintf('End Trial without End stim: %d (%.2f)\n',nt-1,Events.times(j)/10000);
            Trials.End(nt) = Events.times(j);
            Trials.Result(nt) = -1;  % this will be set to 0 if a BadFix is found.
            Trials.TrueEnd(nt) = NaN;
        end
    elseif Events.codes(j,1) == BADFIX & storestate %% Doesn't happen. Badfix is in Text, because it is sent, not received
        Trials.End(nt) = Events.times(j);
        Trials.Result(nt) = 0;
        instim = 0;
    elseif Events.codes(j,1) == STARTEXPT & storestate
        if inexpt %close an existing expt (e.g. if crashed out)
            Expts(nx).end = Events.times(j);
            Expts(nx).lasttrial = nt;
        end
        nx = nx+1;
        Expts(nx).start = Events.times(j);
        Expts(nx).firsttrial = nt+1;
        inexpt = 1;
    elseif Events.codes(j,1) == ENDEXPT & nx & storestate
        Expts(nx).end = Events.times(j);
        Expts(nx).lasttrial = nt;
        inexpt = 0;
        Expts(nx).result = ENDEXPT;
    elseif Events.codes(j,1) == CANCELEXPT & nx & storestate
        Expts(nx).end = Events.times(j);
        Expts(nx).lasttrial = nt;
        Expts(nx).result = CANCELEXPT;
        inexpt = 0;
    elseif Events.codes(j,1) == ENDEXPT
        nx = nx;        
    end
end

Trials.id = zeros(size(Trials.Trial))'; %needs to be a row 
if nt == 0
    Expts = [];
    return;
end
if inexpt
    Expts(nx).lasttrial =nt;
end
ntrials = nt;
fprintf('Setting Trials Took %.2f\n',toc);
