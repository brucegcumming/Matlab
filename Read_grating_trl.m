% Program to read in Bruce's rds.ODX.trl files.
function ORTC = Read_grating_trl(ORTCfile);
fid = fopen(ORTCfile,'r');
if fid==-1
    % file could not be opened
    ORTC=0;
    return
end
ORTC.fileused = ORTCfile;


% ORTuning.data(1) = stim presented to left eye only
% ORTuning.data(2) = stim presented to left eye only
% ORTuning.data(3) = stim presented to both eyes at once
% ORTuning.data(4) = stim presented to neither eye (blank)

ORTC.preperiod.spikecount = [];
ORTC.preperiod.timesincestim = [];
for EYE=1:3
    ORTC.grating(EYE).orientation = [];
    ORTC.grating(EYE).spatfreq = [];
    ORTC.grating(EYE).tempfreq = [];
    ORTC.grating(EYE).disp_horiz = [];
    ORTC.grating(EYE).disp_vert = [];
    ORTC.grating(EYE).disparityphase = [];
    ORTC.grating(EYE).spikecount = [];
    ORTC.grating(EYE).stim_starttime = [];
    ORTC.grating(EYE).stim_duration = [];
    ORTC.grating(EYE).stimwidth = [];
    ORTC.grating(EYE).stimheight = [];
    ORTC.grating(EYE).fitmin = [];
    ORTC.grating(EYE).fitmax = [];
end

while ~feof(fid)
    scrap = fgets(fid);
    
    % Look for line containing word "Preperiod"; assume that first number following this is duration of preperiod in secs
    indx = findstr(lower(scrap),'preperiod');
    if ~isempty(indx)
        ORTC.preperiod.stim_duration = sscanf(scrap(indx+9:end),'%f',1);
    end
    
    % Look for lines beginning with fitmin, fitmax (orientation window to use)
    [word,number] = strtok(scrap);
    if strcmp(word,'Fitmin')
        for EYE=1:3 % fit range is the same for all eyes
            ORTC.grating(EYE).fitmin = str2num(number);
        end
    elseif strcmp(word,'Fitmax')
        for EYE=1:3 % fit range is the same for all eyes
            ORTC.grating(EYE).fitmax = str2num(number);
        end
    end
    
    % Look for lines beginning with a digit
    trialno = sscanf(scrap,'%f');
    if ~isempty(trialno) %then the line begins with a digit
        indx_lbrack = find(scrap=='(');
        indx_rbrack = find(scrap==')');
        indx_comma = find(scrap==',');
        indx_colon = find(scrap==':');
        stim_starttime = 1e4*str2num(scrap((indx_lbrack(1)+1):(indx_comma(1)-1))) ; % convert to seconds
        stim_duration = str2num(scrap((indx_comma(1)+1):(indx_rbrack(1)-1))) ;
        % Select those lines where stim is correlated disparate grating:
        if strcmp(scrap(indx_colon(1)+[2:12]),'grating0-rc')  & ~isempty(findstr(scrap,'cl='))
            stimdescript = scrap((indx_colon(1)+13):end);
            % Remove leading preperiods
            stimdescript=fliplr(deblank(fliplr(stimdescript)));
            % Spike count is first number after the first space in stimdescript:
            [str,tok]=strtok(stimdescript);
            spikecount = sscanf(tok,'%d1'); % save this for later as we do not yet know which eye

            % Spike count during preperiod (preperiod screen) is first number after word "pre":
            indx = findstr(scrap,'pre');
            tmp = sscanf(scrap(indx+3:end),'%f',3); % first entry is spike count in this cluster, 2nd entry is spike count in hash, 3rd entry is time elapsed since last stim
            ORTC.preperiod.spikecount = [ ORTC.preperiod.spikecount tmp(1) ];
            ORTC.preperiod.timesincestim = [ ORTC.preperiod.timesincestim tmp(3) ];
                        
            EYE=[];
            while length(stimdescript)>0
                [param,rem]=strtok(stimdescript,'=');
                val = sscanf(rem(2:end),'%f');
                % Work out what param this value relates to and record it:
                switch strrep(param,',','') %remove any unwanted commas
                case 'lxo'
                    EYE=0;
                case 'rxo'
                    EYE=1;
                case 'bxo'
                    EYE=2;
                case 'or'
                    ORTC.grating(EYE).orientation = [ ORTC.grating(EYE).orientation val];
                case 'sf'
                    ORTC.grating(EYE).spatfreq = [ORTC.grating(EYE).spatfreq val];
                case 'tf'
                    ORTC.grating(EYE).tempfreq = [ORTC.grating(EYE).tempfreq val];
                case 'dx'
                    ORTC.grating(EYE).disp_horiz = [ ORTC.grating(EYE).disp_horiz val];
                case 'dy'
                    ORTC.grating(EYE).disp_vert = [ORTC.grating(EYE).disp_vert val];
                case 'dp'
                    ORTC.grating(EYE).disparityphase = [ORTC.grating(EYE).disparityphase val];
                case 'wi'
                    ORTC.grating(EYE).stimwidth =  [ORTC.grating(EYE).stimwidth val];
                case 'hi'
                    ORTC.grating(EYE).stimheight =  [ORTC.grating(EYE).stimheight val];
                end
                % Go up to next comma:
                [rem,stimdescript] = strtok(rem,',');
            end % finish going through stim description
            % Record the spikecount
            ORTC.grating(EYE).spikecount = [ORTC.grating(EYE).spikecount spikecount];
            ORTC.grating(EYE).stim_starttime = [ORTC.grating(EYE).stim_starttime stim_starttime];
            ORTC.grating(EYE).stim_duration = [ORTC.grating(EYE).stim_duration stim_duration];
        end % if grating
                
        % Select those lines where stim is blanl:
        if strcmp(scrap(indx_colon(1)+[2:12]),'none0-rc')
            stimdescript = scrap((indx_colon(1)+9):end);
            % Remove leading blank
            stimdescript=fliplr(deblank(fliplr(stimdescript)));
            % Spike count during stimulus is first number after the first space in stimdescript:
            [str,tok]=strtok(stimdescript);
            spikecount = sscanf(tok,'%d1'); % save this for later as we do not yet know which eye

            % Spike count during preperiod (preperiod screen) is first number after word "pre":
            indx = findstr(scrap,'pre');
            tmp = sscanf(scrap(indx+3:end),'%f',3); % first entry is spike count in this cluster, 2nd entry is spike count in hash, 3rd entry is time elapsed since last stim
            ORTC.preperiod.spikecount = [ ORTC.preperiod.spikecount tmp(1) ];
            ORTC.preperiod.timesincestim = [ ORTC.preperiod.timesincestim tmp(3) ];
            
            % Record the spikecount
            ORTC.blank.spikecount = [ORTC.blank.spikecount spikecount];
            ORTC.blank.stim_starttime = [ORTC.blank.stim_starttime stim_starttime];
            ORTC.blank.stim_duration = [ORTC.blank.stim_duration stim_duration];
        end % if blank
        
    end
end

% Throw away any data where the stim duration was 10% > than average, as this is probably a cock-up
for EYE=1:length(ORTC.grating)
    if isfield(ORTC.grating(EYE),'stim_duration')
        mn = mean(ORTC.grating(EYE).stim_duration);
        nn = length(ORTC.grating(EYE).stim_duration);
        indx = find(0 < ORTC.grating(EYE).stim_duration & ORTC.grating(EYE).stim_duration < 1.1*mn);
        % Go through, look for all fields which have the same length as the stim_duration field,
        % and select only those which satisfy the criteria
        agri = fieldnames(ORTC.grating(EYE));
        nagri = length(agri);
        for jager = 1:nagri
            ager = getfield(ORTC.grating(EYE),agri{jager});
            if length(ager)==nn
                ORTC.grating(EYE) = setfield(ORTC.grating(EYE),agri{jager},ager(indx));
            end
        end
    end
end

% For the preperiod, only use data recorded in sessions where the time since last stim was between 0s and some maximum, specified in the DEFINITIONS file
if ~isempty(ORTC.preperiod.spikecount)
    ORTC.preperiod.firingrate = ORTC.preperiod.spikecount ./ ORTC.preperiod.stim_duration;
    ORTC.preperiod.avgd.indxtoav = find(ORTC.preperiod.timesincestim > MINTIMEELAPSED & ORTC.preperiod.timesincestim < MAXTIMEELAPSED );
    ORTC.preperiod.avgd.firingrates = ORTC.preperiod.firingrate(ORTC.preperiod.avgd.indxtoav)
    ORTC.preperiod.avgd.firingrate_mean = mean(ORTC.preperiod.avgd.firingrates);
    ORTC.preperiod.avgd.firingrate_SD = std(ORTC.preperiod.avgd.firingrates);
    ORTC.preperiod.avgd.nrep = length(ORTC.preperiod.avgd.indxtoav);
    ORTC.preperiod.avgd.firingrate_SEM = ORTC.preperiod.avgd.firingrate_SD./sqrt(ORTC.preperiod.avgd.nrep);
end

% Now go through and take the averages
for EYE=1:3
    ORTC.grating(EYE).firingrate = ORTC.grating(EYE).spikecount ./ ORTC.grating(EYE).stim_duration; % spikes per sec    
    ORTC.grating(EYE).nOR = length(ORTC.grating(EYE).orientation);    
    
    if ORTC.grating(EYE).nOR>1
        % Work out min and max
        if isempty(ORTC.grating(EYE).fitmin)
            ORTC.grating(EYE).fitmin = min(ORTC.grating(EYE).orientation);
        end
        if isempty(ORTC.grating(EYE).fitmax)
            ORTC.grating(EYE).fitmax = max(ORTC.grating(EYE).orientation);
        end
        % Alert the user if fitrange!=180:
        if ORTC.grating(EYE).fitmax-ORTC.grating(EYE).fitmin ~= 180
            ORTC.grating(EYE).comment = 'fitting range is not 180 degrees ... something wrong?';
        end
        
        % First just average everything that was recorded at the same spatfreq, whether or not other params were identical
        distinctOR = Members(ORTC.grating(EYE).spatfreq);
        distinctOR = sort(distinctOR);
        for j=1:length(distinctOR)
            indx = find( ORTC.grating(EYE).spatfreq==distinctOR(j));
            ORTC.grating(EYE).allavgd.firingrate_mean(j) = mean(ORTC.grating(EYE).firingrate(indx));
            ORTC.grating(EYE).allavgd.firingrates{j} = ORTC.grating(EYE).firingrate(indx);
            ORTC.grating(EYE).allavgd.firingrate_SD(j) = std(ORTC.grating(EYE).firingrate(indx));
            ORTC.grating(EYE).allavgd.nrep(j) = length(indx);
            ORTC.grating(EYE).allavgd.firingrate_meansqrt(j) = mean(sqrt(ORTC.grating(EYE).firingrate(indx)));
        end
        ORTC.grating(EYE).allavgd.firingrate_sqmeansqrt = (ORTC.grating(EYE).allavgd.firingrate_meansqrt).^2;
        ORTC.grating(EYE).allavgd.firingrate_SEM = ORTC.grating(EYE).allavgd.firingrate_SD ./ sqrt(ORTC.grating(EYE).allavgd.nrep);
        ORTC.grating(EYE).allavgd.spatfreq = distinctOR;
        % Work out area under ORTC: Just do a trapezoidal integration 
        if length(distinctOR)>1
            ORTC.grating(EYE).allavgd.area = trapz(distinctOR,ORTC.grating(EYE).allavgd.firingrate_mean);
        end   
        
        % Now Go through it and take the mean of all trials in which stim params were identical
        trialstolookat = [1:ORTC.grating(EYE).nOR];
        nORTCs=0; % count the number of distinct ORTCs contained in this file (Eg, there may be two with different spatfreqs)
        while ~isempty(trialstolookat)
            jj=trialstolookat(1); % compare everything to this
            indx_thisTC = trialstolookat;
            % NB comment the above line and uncomment what follows to make the program check that everything is consistent, and output several TCs if not.
            %indx_thisTC = find( ORTC.grating(EYE).spatfreq==ORTC.grating(EYE).spatfreq(jj) & ORTC.grating(EYE).tempfreq==ORTC.grating(EYE).tempfreq(jj) & ORTC.grating(EYE).disp_horiz==ORTC.grating(EYE).disp_horiz(jj) ...
            %    & ORTC.grating(EYE).disp_vert==ORTC.grating(EYE).disp_vert(jj)  & ORTC.grating(EYE).disparityphase==ORTC.grating(EYE).disparityphase(jj) ...
            %    & ORTC.grating(EYE).stimwidth==ORTC.grating(EYE).stimwidth(jj)  & ORTC.grating(EYE).stimheight==ORTC.grating(EYE).stimheight(jj) );
            if ~isempty(indx_thisTC)
                nORTCs = nORTCs+1;
                % Work out what distinct ORs there are in this TC
                % Extract distinct ORs: record num rep, mean and SD
                distinctOR = Members(ORTC.grating(EYE).orientation(indx_thisTC));
                distinctOR = sort(distinctOR);
                % Record the params used in getting this TC:
                ORTC.grating(EYE).avgd(nORTCs).spatfreq = ORTC.grating(EYE).spatfreq(jj);
                ORTC.grating(EYE).avgd(nORTCs).tempfreq = ORTC.grating(EYE).tempfreq(jj);
                ORTC.grating(EYE).avgd(nORTCs).disp_horiz = ORTC.grating(EYE).disp_horiz(jj);
                ORTC.grating(EYE).avgd(nORTCs).disp_vert = ORTC.grating(EYE).disp_vert(jj);
                ORTC.grating(EYE).avgd(nORTCs).disparityphase = ORTC.grating(EYE).disparityphase(jj);
                ORTC.grating(EYE).avgd(nORTCs).stimwidth = ORTC.grating(EYE).stimwidth(jj);
                ORTC.grating(EYE).avgd(nORTCs).stimheight = ORTC.grating(EYE).stimheight(jj);
                
                % take the mean of everything:
                for j=1:length(distinctOR)
                    indx = find( ORTC.grating(EYE).orientation==distinctOR(j) & ORTC.grating(EYE).spatfreq==ORTC.grating(EYE).spatfreq(jj) & ORTC.grating(EYE).tempfreq==ORTC.grating(EYE).tempfreq(jj) ...
                        & ORTC.grating(EYE).disp_horiz==ORTC.grating(EYE).disp_horiz(jj)  & ORTC.grating(EYE).disp_vert==ORTC.grating(EYE).disp_vert(jj)  & ORTC.grating(EYE).disparityphase==ORTC.grating(EYE).disparityphase(jj) ...
                        & ORTC.grating(EYE).stimwidth==ORTC.grating(EYE).stimwidth(jj)  & ORTC.grating(EYE).stimheight==ORTC.grating(EYE).stimheight(jj) );
                    ORTC.grating(EYE).avgd(nORTCs).firingrates{j} = ORTC.grating(EYE).firingrate(indx);
                    ORTC.grating(EYE).avgd(nORTCs).firingrate_mean(j) = mean(ORTC.grating(EYE).firingrate(indx));
                    ORTC.grating(EYE).avgd(nORTCs).firingrate_SD(j) = std(ORTC.grating(EYE).firingrate(indx));
                    ORTC.grating(EYE).avgd(nORTCs).nrep(j) = length(indx);
                    ORTC.grating(EYE).avgd(nORTCs).firingrate_meansqrt(j) = mean(sqrt(ORTC.grating(EYE).firingrate(indx)));
                end
                ORTC.grating(EYE).avgd(nORTCs).firingrate_sqmeansqrt = (ORTC.grating(EYE).avgd(nORTCs).firingrate_meansqrt).^2;
                
                ORTC.grating(EYE).avgd(nORTCs).firingrate_SEM = ORTC.grating(EYE).avgd(nORTCs).firingrate_SD ./ sqrt(ORTC.grating(EYE).avgd(nORTCs).nrep);
                ORTC.grating(EYE).avgd(nORTCs).orientation = distinctOR;
                
                % Work out area under OTC: Just do a trapezoidal integration between fitmin and fitmax
                if length(distinctOR)>1
                    indx = find(distinctOR>=ORTC.grating(EYE).fitmin & distinctOR<=ORTC.grating(EYE).fitmax);
                    ORTC.grating(EYE).avgd(nORTCs).area = trapz(distinctOR(indx),ORTC.grating(EYE).avgd(nORTCs).firingrate_mean(indx));
                end
                
            end
            % remove those trials which have now been averaged into some TC:
            for j=1:length(indx_thisTC)
                trialstolookat = trialstolookat(find(trialstolookat~=indx_thisTC(j)));
            end
            
        end % next tuning curve
    end % checking that there is any data for this eye
end % next eye

fclose(fid)