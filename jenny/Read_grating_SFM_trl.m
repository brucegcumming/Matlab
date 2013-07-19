% Program to read in Bruce's rds.Odisp_horiz.trl files.
function SFTC = Read_grating_SFM_Trl(SFTCfile);
DEFINITIONS
fid = fopen(SFTCfile,'r');
if fid==-1
    % file could not be opened
    SFTC=0;
    return
end
SFTC.fileused = SFTCfile;

SFTC.preperiod.spikecount = [];
SFTC.preperiod.timesincestim = [];
for EYE=1:3
    SFTC.grating(EYE).orientation = [];
    SFTC.grating(EYE).spatfreq = [];
    SFTC.grating(EYE).tempfreq = [];
    SFTC.grating(EYE).disp_horiz = [];
    SFTC.grating(EYE).disp_vert = [];
    SFTC.grating(EYE).disparityphase = [];
    SFTC.grating(EYE).spikecount = [];
    SFTC.grating(EYE).stim_starttime = [];
    SFTC.grating(EYE).stim_duration = [];
    SFTC.grating(EYE).stimwidth = [];
    SFTC.grating(EYE).stimheight = [];
end

while ~feof(fid)
    scrap = fgets(fid);
    
    % Look for line containing word "preperiod"; assume that first number following this is duration of preperiod in secs
    indx = findstr(lower(scrap),'preperiod');
    if ~isempty(indx)
        SFTC.preperiod.stim_duration = sscanf(scrap(indx+9:end),'%f',1);
    end
    
    % Look for lines beginning with digit (trial)
    trialno = sscanf(scrap,'%f');
    if ~isempty(trialno) %then the line begins with a digit
        indx_lbrack = find(scrap=='(');
        indx_rbrack = find(scrap==')');
        indx_comma = find(scrap==',');
        indx_colon = find(scrap==':');
        stim_starttime = 1e4*str2num(scrap((indx_lbrack(1)+1):(indx_comma(1)-1))) ; % convert to seconds
        stim_duration = str2num(scrap((indx_comma(1)+1):(indx_rbrack(1)-1))) ;
        
        % Select those lines where stim is correlated disparate grating:
        if strcmp(scrap(indx_colon(1)+[2:12]),'grating0-rc') & ~isempty(findstr(scrap,'cl='))
            stimdescript = scrap((indx_colon(1)+13):end);
            % Remove leading blank
            stimdescript=fliplr(deblank(fliplr(stimdescript)));
            % Spike count during stimulus is first number after the first space in stimdescript:
            [str,tok]=strtok(stimdescript);
            spikecount = sscanf(tok,'%d1'); % save this for later as we do not yet know which eye
            
            % Spike count during preperiod (preperiod screen) is first number after word "pre":
            indx = findstr(scrap,'pre');
            tmp = sscanf(scrap(indx+3:end),'%f',3); % first entry is spike count in this cluster, 2nd entry is spike count in hash, 3rd entry is time elapsed since last stim
            SFTC.preperiod.spikecount = [ SFTC.preperiod.spikecount tmp(1) ];
            SFTC.preperiod.timesincestim = [ SFTC.preperiod.timesincestim tmp(3) ];
            
            EYE=[];
            while length(stimdescript)>0
                [param,rem]=strtok(stimdescript,'=');
                val = sscanf(rem(2:end),'%f');
                % Work out what param this value relates to and record it:
                switch strrep(param,',','') %remove any unwanted commas
                case 'lxo'
                    EYE=L;
                case 'rxo'
                    EYE=R;
                case 'bxo'
                    EYE=B;
                case 'or'
                    SFTC.grating(EYE).orientation = [ SFTC.grating(EYE).orientation val];
                case 'sf'
                    SFTC.grating(EYE).spatfreq = [SFTC.grating(EYE).spatfreq val];
                case 'tf'
                    SFTC.grating(EYE).tempfreq = [SFTC.grating(EYE).tempfreq val];
                case 'dx'
                    SFTC.grating(EYE).disp_horiz = [ SFTC.grating(EYE).disp_horiz val];
                case 'dy'
                    SFTC.grating(EYE).disp_vert = [SFTC.grating(EYE).disp_vert val];
                case 'dp'
                    SFTC.grating(EYE).disparityphase = [SFTC.grating(EYE).disparityphase val];
                case 'wi'
                    SFTC.grating(EYE).stimwidth =  [SFTC.grating(EYE).stimwidth val];
                case 'hi'
                    SFTC.grating(EYE).stimheight =  [SFTC.grating(EYE).stimheight val];
                end
                % Go up to next comma:
                [rem,stimdescript] = strtok(rem,',');
            end % finish going through stim description
            % Record the spikecount
            SFTC.grating(EYE).spikecount = [SFTC.grating(EYE).spikecount spikecount];
            SFTC.grating(EYE).stim_starttime = [SFTC.grating(EYE).stim_starttime stim_starttime];
            SFTC.grating(EYE).stim_duration = [SFTC.grating(EYE).stim_duration stim_duration];
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
            SFTC.preperiod.spikecount = [ SFTC.preperiod.spikecount tmp(1) ];
            SFTC.preperiod.timesincestim = [ SFTC.preperiod.timesincestim tmp(3) ];
            
            % Record the spikecount
            SFTC.blank.spikecount = [SFTC.blank.spikecount spikecount];
            SFTC.blank.stim_starttime = [SFTC.blank.stim_starttime stim_starttime];
            SFTC.blank.stim_duration = [SFTC.blank.stim_duration stim_duration];
        end % if blank        
    end
end

% Throw away any data where the stim duration was 10% > than average, as this is probably a cock-up.
% Also throw away any data recorded at spat freq of > 20cpd - this is unreliable as it will probably alias
for EYE=1:length(SFTC.grating)
    if isfield(SFTC.grating(EYE),'stim_duration')
        if ~isempty(SFTC.grating(EYE).stim_duration)
            mn = mean(SFTC.grating(EYE).stim_duration);
            nn = length(SFTC.grating(EYE).stim_duration);
            indx = find(0 < SFTC.grating(EYE).stim_duration & SFTC.grating(EYE).stim_duration < 1.1*mn & SFTC.grating(EYE).spatfreq <= MAXALLOWEDSF);
            % Go through, look for all fields which have the same length as the stim_duration field,
            % and select only those which satisfy the criteria
            agri = fieldnames(SFTC.grating(EYE));
            nagri = length(agri);
            for jager = 1:nagri
                ager = getfield(SFTC.grating(EYE),agri{jager});
                if length(ager)==nn
                    SFTC.grating(EYE) = setfield(SFTC.grating(EYE),agri{jager},ager(indx));
                end
            end
        end
    end
end

if ~isempty(SFTC.preperiod.spikecount)
    SFTC.preperiod.firingrate = SFTC.preperiod.spikecount ./ SFTC.preperiod.stim_duration;
    SFTC.preperiod.avgd.indxtoav = find(SFTC.preperiod.timesincestim > MINTIMEELAPSED & SFTC.preperiod.timesincestim < MAXTIMEELAPSED );
    SFTC.preperiod.avgd.firingrates = SFTC.preperiod.firingrate(SFTC.preperiod.avgd.indxtoav);
    SFTC.preperiod.avgd.firingrate_mean = mean(SFTC.preperiod.avgd.firingrates);
    SFTC.preperiod.avgd.firingrate_SD = std(SFTC.preperiod.avgd.firingrates);
    SFTC.preperiod.avgd.nrep = length(SFTC.preperiod.avgd.indxtoav);
    SFTC.preperiod.avgd.firingrate_SEM = SFTC.preperiod.avgd.firingrate_SD./sqrt(SFTC.preperiod.avgd.nrep);
end

% Average to construct actual tuning curves
for EYE=1:3
    SFTC.grating(EYE).firingrate = SFTC.grating(EYE).spikecount ./ SFTC.grating(EYE).stim_duration; % spikes per sec    
    SFTC.grating(EYE).nSF = length(SFTC.grating(EYE).spatfreq);
    
    if SFTC.grating(EYE).nSF>0        
        
        % First just average everything that was recorded at the same spatfreq, whether or not other params were identical
        distinctSF = Members(SFTC.grating(EYE).spatfreq);
        distinctSF = sort(distinctSF);
        for j=1:length(distinctSF)
            indx = find( SFTC.grating(EYE).spatfreq==distinctSF(j));
            SFTC.grating(EYE).allavgd.firingrate_mean(j) = mean(SFTC.grating(EYE).firingrate(indx));
            SFTC.grating(EYE).allavgd.firingrates{j} = SFTC.grating(EYE).firingrate(indx);
            SFTC.grating(EYE).allavgd.firingrate_SD(j) = std(SFTC.grating(EYE).firingrate(indx));
            SFTC.grating(EYE).allavgd.nrep(j) = length(indx);
            SFTC.grating(EYE).allavgd.firingrate_meansqrt(j) = mean(sqrt(SFTC.grating(EYE).firingrate(indx)));
        end
        SFTC.grating(EYE).allavgd.firingrate_sqmeansqrt = (SFTC.grating(EYE).allavgd.firingrate_meansqrt).^2;
        SFTC.grating(EYE).allavgd.firingrate_SEM = SFTC.grating(EYE).allavgd.firingrate_SD ./ sqrt(SFTC.grating(EYE).allavgd.nrep);
        SFTC.grating(EYE).allavgd.spatfreq = distinctSF;
        % Work out area under SFTC: Just do a trapezoidal integration 
        if length(distinctSF)>1
            SFTC.grating(EYE).allavgd.area = trapz(distinctSF,SFTC.grating(EYE).allavgd.firingrate_mean);
        end        
        
        % Now Go through it and take the mean of all trials in which stim params were identical
        trialstolookat = [1:SFTC.grating(EYE).nSF];
        nSFTCs=0; % count the number of distinct SFTCs contained in this file (Eg, there may be two with different orientations)
        while ~isempty(trialstolookat)
            jj=trialstolookat(1); % compare everything to this
            indx_thisTC = find( SFTC.grating(EYE).orientation==SFTC.grating(EYE).orientation(jj) & SFTC.grating(EYE).tempfreq==SFTC.grating(EYE).tempfreq(jj) & SFTC.grating(EYE).disp_horiz==SFTC.grating(EYE).disp_horiz(jj) ...
                & SFTC.grating(EYE).disp_vert==SFTC.grating(EYE).disp_vert(jj)  & SFTC.grating(EYE).disparityphase==SFTC.grating(EYE).disparityphase(jj) ...
                & SFTC.grating(EYE).stimwidth==SFTC.grating(EYE).stimwidth(jj)  & SFTC.grating(EYE).stimheight==SFTC.grating(EYE).stimheight(jj) );
            if ~isempty(indx_thisTC)
                nSFTCs = nSFTCs+1;
                % Work out what distinct SFs there are in this TC
                % Extract distinct SFs: record num rep, mean and SD
                distinctSF = Members(SFTC.grating(EYE).spatfreq(indx_thisTC));
                distinctSF = sort(distinctSF);
                % Record the params used in getting this TC:
                SFTC.grating(EYE).avgd(nSFTCs).orientation = SFTC.grating(EYE).orientation(jj);
                SFTC.grating(EYE).avgd(nSFTCs).tempfreq = SFTC.grating(EYE).tempfreq(jj);
                SFTC.grating(EYE).avgd(nSFTCs).disparityphase = SFTC.grating(EYE).disparityphase(jj);
                SFTC.grating(EYE).avgd(nSFTCs).disp_horiz = SFTC.grating(EYE).disp_horiz(jj);
                SFTC.grating(EYE).avgd(nSFTCs).disp_vert = SFTC.grating(EYE).disp_vert(jj);
                SFTC.grating(EYE).avgd(nSFTCs).stimwidth = SFTC.grating(EYE).stimwidth(jj);
                SFTC.grating(EYE).avgd(nSFTCs).stimheight = SFTC.grating(EYE).stimheight(jj);
                
                % take the mean of everything:
                for j=1:length(distinctSF)
                    indx = find( SFTC.grating(EYE).spatfreq==distinctSF(j) & SFTC.grating(EYE).orientation==SFTC.grating(EYE).orientation(jj) & SFTC.grating(EYE).tempfreq==SFTC.grating(EYE).tempfreq(jj) ...
                        & SFTC.grating(EYE).disp_horiz==SFTC.grating(EYE).disp_horiz(jj)  & SFTC.grating(EYE).disp_vert==SFTC.grating(EYE).disp_vert(jj)  & SFTC.grating(EYE).disparityphase==SFTC.grating(EYE).disparityphase(jj) ...
                        & SFTC.grating(EYE).stimwidth==SFTC.grating(EYE).stimwidth(jj)  & SFTC.grating(EYE).stimheight==SFTC.grating(EYE).stimheight(jj) );
                    SFTC.grating(EYE).avgd(nSFTCs).firingrate_mean(j) = mean(SFTC.grating(EYE).firingrate(indx));
                    SFTC.grating(EYE).avgd(nSFTCs).firingrates{j} = SFTC.grating(EYE).firingrate(indx);
                    SFTC.grating(EYE).avgd(nSFTCs).firingrate_SD(j) = std(SFTC.grating(EYE).firingrate(indx));
                    SFTC.grating(EYE).avgd(nSFTCs).nrep(j) = length(indx);
                    SFTC.grating(EYE).avgd(nSFTCs).firingrate_meansqrt(j) = mean(sqrt(SFTC.grating(EYE).firingrate(indx)));
                end
                SFTC.grating(EYE).avgd(nSFTCs).firingrate_sqmeansqrt = (SFTC.grating(EYE).avgd(nSFTCs).firingrate_meansqrt).^2;
                SFTC.grating(EYE).avgd(nSFTCs).firingrate_SEM = SFTC.grating(EYE).avgd(nSFTCs).firingrate_SD ./ sqrt(SFTC.grating(EYE).avgd(nSFTCs).nrep);
                SFTC.grating(EYE).avgd(nSFTCs).spatfreq = distinctSF;
                % Work out area under SFTC: Just do a trapezoidal integration 
                if length(distinctSF)>1
                    SFTC.grating(EYE).avgd(nSFTCs).area = trapz(distinctSF,SFTC.grating(EYE).avgd(nSFTCs).firingrate_mean);
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
