% Program to read in Bruce's rds.ODX.trl files.
function RDSresponse = Read_rds_ODX_trl(DTCfile);
DEFINITIONS
fid = fopen(DTCfile,'r');
if fid==-1
    % file could not be opened
    RDSresponse=0;
    return
end
DTC.fileused = DTCfile;

DTC.disp_horiz = [];
DTC.disp_vert = [];
DTC.disp_orthog = [];
DTC.disp_paral = [];
DTC.dotwidth = [];
DTC.spikecount = [];
DTC.stimwidth = [];
DTC.stimheight = [];
DTC.stim_starttime = [];
DTC.stim_duration = [];
DTC.comment=[];

uncorr.stim_starttime = [];
uncorr.stim_duration = [];
uncorr.dotwidth = [];
uncorr.spikecount = [];
uncorr.stimwidth = [];
uncorr.stimheight = [];
uncorr.comment=[];

Lmonoc.stim_starttime = [];
Lmonoc.stim_duration = [];
Lmonoc.dotwidth = [];
Lmonoc.spikecount = [];
Lmonoc.stimwidth = [];
Lmonoc.stimheight = [];
Lmonoc.comment=[];

Rmonoc.stim_starttime = [];
Rmonoc.stim_duration = [];
Rmonoc.dotwidth = [];
Rmonoc.spikecount = [];
Rmonoc.stimwidth = [];
Rmonoc.stimheight = [];
Rmonoc.comment=[];

blank.stim_starttime = [];
blank.stim_duration = [];
blank.spikecount = [];
blank.comment=[];

preperiod.spikecount=[];
preperiod.timesincestim=[];

while ~feof(fid)
    scrap = fgets(fid);
    
    % Look for line containing word "Preperiod"; assume that first number following this is duration of preperiod in secs
    indx = findstr(lower(scrap),'preperiod');
    if ~isempty(indx)
        preperiod.stim_duration = sscanf(scrap(indx+9:end),'%f',1);
    end
    
    trialno = sscanf(scrap,'%f');
    if ~isempty(trialno) %then the line begins with a digit
        indx_lbrack = find(scrap=='(');
        indx_rbrack = find(scrap==')');
        indx_comma = find(scrap==',');
        indx_colon = find(scrap==':');
        stim_starttime = 1e4*str2num(scrap((indx_lbrack(1)+1):(indx_comma(1)-1))) ; % convert to seconds
        stim_duration = str2num(scrap((indx_comma(1)+1):(indx_rbrack(1)-1))) ;
        
        % Spike count during preperiod (preperiod screen) is first number after word "pre":
        indx = findstr(scrap,'pre');
        tmp = sscanf(scrap(indx+3:end),'%f',3); % first entry is spike count in this cluster, 2nd entry is spike count in hash, 3rd entry is time elapsed since last stim
        preperiod.spikecount = [ preperiod.spikecount tmp(1) ];
        preperiod.timesincestim = [ preperiod.timesincestim tmp(3) ];
        
        %%%%%%%%%%%%%%%%%% stim is correlated disparate RDS:
        if strcmp(scrap(indx_colon(1)+[2:8]),'rds0-rc') & isempty(findstr(scrap,'+uc')) & ~isempty(findstr(scrap,'bxo')) & ~isempty(findstr(scrap,'cl='))
            DTC.stim_starttime = [DTC.stim_starttime stim_starttime];
            DTC.stim_duration = [DTC.stim_duration stim_duration];
            stimdescript = scrap((indx_colon(1)+9):end);
            % Remove leading blanks
            stimdescript=fliplr(deblank(fliplr(stimdescript)));
            % Spike count is first number after the first space in stimdescript:
            [str,tok]=strtok(stimdescript);
            DTC.spikecount = [DTC.spikecount sscanf(tok,'%d1')];
            
            while length(stimdescript)>0
                [param,rem]=strtok(stimdescript,'=');
                val = sscanf(rem(2:end),'%f');
                % Work out what param this value relates to and record it:
                switch strrep(param,',','') %remove any unwanted commas
                case 'dx'
                    DTC.disp_horiz = [ DTC.disp_horiz val];
                case 'dy'
                    DTC.disp_vert = [ DTC.disp_vert val];
                case 'dO'
                    DTC.disp_orthog = [ DTC.disp_orthog val];
                case 'dP'
                    DTC.disp_paral = [ DTC.disp_paral val];
                case 'dw'
                    DTC.dotwidth = [DTC.dotwidth val];
                case 'wi'
                    DTC.stimwidth =  [DTC.stimwidth val];
                case 'hi'
                    DTC.stimheight =  [DTC.stimheight val];
                end
                % Go up to next comma:
                [rem,stimdescript] = strtok(rem,',');
            end % finish going through stim description
            
        end % if correlated disparate RDS
        
        %%%%%%%%%%%%%%%%%% stim is binocularly uncorrelated disparate RDS:
        if strcmp(scrap(indx_colon(1)+[2:8]),'rds0-rc') & ~isempty(findstr(scrap,'+uc')) & ~isempty(findstr(scrap,'bxo')) & ~isempty(findstr(scrap,'cl='))
            uncorr.stim_starttime = [uncorr.stim_starttime stim_starttime];
            uncorr.stim_duration = [uncorr.stim_duration stim_duration];
            stimdescript = scrap((indx_colon(1)+9):end);
            % Remove leading blanks
            stimdescript=fliplr(deblank(fliplr(stimdescript)));
            % Spike count is first number after the first space in stimdescript:
            [str,tok]=strtok(stimdescript);
            uncorr.spikecount = [uncorr.spikecount sscanf(tok,'%d1')];
            while length(stimdescript)>0
                [param,rem]=strtok(stimdescript,'=');
                val = sscanf(rem(2:end),'%f');
                % Work out what param this value relates to and record it:
                switch strrep(param,',','') %remove any unwanted commas
                case 'dw'
                    uncorr.dotwidth = [uncorr.dotwidth val];
                case 'wi'
                    uncorr.stimwidth =  [uncorr.stimwidth val];
                case 'hi'
                    uncorr.stimheight =  [uncorr.stimheight val];
                end
                % Go up to next comma:
                [rem,stimdescript] = strtok(rem,',');
            end % finish going through stim description
            
        end % if uncorrelated RDS
        
        %%%%%%%%%%%%%%%%%% stim is blank:
        if strcmp(scrap(indx_colon(1)+[2:9]),'none0-rc')  & ~isempty(findstr(scrap,'cl='))
            blank.stim_starttime = [blank.stim_starttime stim_starttime];
            blank.stim_duration = [blank.stim_duration stim_duration];
            stimdescript = scrap((indx_colon(1)+10):end);
            % Remove leading blanks
            stimdescript=fliplr(deblank(fliplr(stimdescript)));
            % Spike count is first number after the first space in stimdescript:
            [str,tok]=strtok(stimdescript);
            blank.spikecount = [blank.spikecount sscanf(tok,'%d1')];            
        end % if blank
        
        %%%%%%%%%%%%%%%%%% stim is left monocular RDS:
        if strcmp(scrap(indx_colon(1)+[2:8]),'rds0-rc') & ~isempty(findstr(scrap,'lxo'))  & ~isempty(findstr(scrap,'cl='))
            Lmonoc.stim_starttime = [Lmonoc.stim_starttime stim_starttime];
            Lmonoc.stim_duration = [Lmonoc.stim_duration stim_duration];
            stimdescript = scrap((indx_colon(1)+9):end);
            % Remove leading blanks
            stimdescript=fliplr(deblank(fliplr(stimdescript)));
            % Spike count is first number after the first space in stimdescript:
            [str,tok]=strtok(stimdescript);
            Lmonoc.spikecount = [Lmonoc.spikecount sscanf(tok,'%d1')];
            while length(stimdescript)>0
                [param,rem]=strtok(stimdescript,'=');
                val = sscanf(rem(2:end),'%f');
                % Work out what param this value relates to and record it:
                switch strrep(param,',','') %remove any unwanted commas
                case 'dw'
                    Lmonoc.dotwidth = [Lmonoc.dotwidth val];
                case 'wi'
                    Lmonoc.stimwidth =  [Lmonoc.stimwidth val];
                case 'hi'
                    Lmonoc.stimheight =  [Lmonoc.stimheight val];
                end
                % Go up to next comma:
                [rem,stimdescript] = strtok(rem,',');
            end % finish going through stim description            
        end % if L monoc RDS
        
        %%%%%%%%%%%%%%%%%% stim is right monocular RDS:
        if strcmp(scrap(indx_colon(1)+[2:8]),'rds0-rc') & ~isempty(findstr(scrap,'rxo'))  & ~isempty(findstr(scrap,'cl='))
            Rmonoc.stim_starttime = [Rmonoc.stim_starttime stim_starttime];
            Rmonoc.stim_duration = [Rmonoc.stim_duration stim_duration];
            stimdescript = scrap((indx_colon(1)+9):end);
            % Remove leading blanks
            stimdescript=fliplr(deblank(fliplr(stimdescript)));
            % Spike count is first number after the first space in stimdescript:
            [str,tok]=strtok(stimdescript);
            Rmonoc.spikecount = [Rmonoc.spikecount sscanf(tok,'%d1')];
            while length(stimdescript)>0
                [param,rem]=strtok(stimdescript,'=');
                val = sscanf(rem(2:end),'%f');
                % Work out what param this value relates to and record it:
                switch strrep(param,',','') %remove any unwanted commas
                case 'dw'
                    Rmonoc.dotwidth = [Rmonoc.dotwidth val];
                case 'wi'
                    Rmonoc.stimwidth =  [Rmonoc.stimwidth val];
                case 'hi'
                    Rmonoc.stimheight =  [Rmonoc.stimheight val];
                end
                % Go up to next comma:
                [rem,stimdescript] = strtok(rem,',');
            end % finish going through stim description            
        end % if R monoc RDS
        
    end
end
fclose(fid);

% Throw away any data where the stim duration was 10% > than average, as this is probably a cock-up
for stimtypecell={'blank' 'Lmonoc' 'Rmonoc' 'DTC' 'uncorr' }
    expt = eval(stimtypecell{1});
    if isfield(expt,'stim_duration')        
        stimdur = getfield(expt,'stim_duration');
        mn = mean(stimdur);
        nn = length(stimdur);
        indx = find(0 < stimdur & stimdur < 1.1*mn);
        if length(indx) < length(stimdur)
            % Go through, look for all fields which have the same length as the stim_duration field,
            % and select only those which satisfy the criteria
            agri = fieldnames(expt);
            nagri = length(agri);
            for jager = 1:nagri
                ager = getfield(expt,agri{jager});
                if length(ager)==nn
                    expt = setfield(expt,agri{jager},ager(indx));                
                end
            end
        end
    end
    eval([stimtypecell{1} '=expt;']);
end

DTC.firingrate = DTC.spikecount ./ DTC.stim_duration; % spikes per sec
Lmonoc.firingrate = Lmonoc.spikecount ./ Lmonoc.stim_duration; % spikes per sec
Rmonoc.firingrate = Rmonoc.spikecount ./ Rmonoc.stim_duration; % spikes per sec
uncorr.firingrate = uncorr.spikecount ./ uncorr.stim_duration; % spikes per sec
blank.firingrate = blank.spikecount ./ blank.stim_duration; % spikes per sec

preperiod.firingrate = preperiod.spikecount ./ preperiod.stim_duration;
preperiod.avgd.indxtoav = find(preperiod.timesincestim > MINTIMEELAPSED & preperiod.timesincestim < MAXTIMEELAPSED);
preperiod.avgd.firingrates = preperiod.firingrate(preperiod.avgd.indxtoav);
preperiod.avgd.firingrate_mean = mean(preperiod.avgd.firingrates);
preperiod.avgd.firingrate_SD = std(preperiod.avgd.firingrates);
preperiod.avgd.nrep = length(preperiod.avgd.indxtoav);
preperiod.avgd.firingrate_SEM = preperiod.avgd.firingrate_SD./sqrt(preperiod.avgd.nrep);

% For L,R,blank and uncorr, average the firing rates to get the mean response.
if ~isempty(blank.spikecount)
    blank.avgd.firingrate_mean = mean(blank.firingrate);
    blank.avgd.firingrate_SD = std(blank.firingrate);
    blank.avgd.firingrate_sqmeansqrt = mean(sqrt(blank.firingrate))^2;
    blank.avgd.nrep = length(blank.firingrate);
    blank.avgd.firingrate_SEM = blank.avgd.firingrate_SD./sqrt(blank.avgd.nrep);
end
% For L,R, uncorr, Record if stim params varied
for stimtypecell = {'Lmonoc' 'Rmonoc' 'uncorr'};
    stimtype = stimtypecell{1};
    dw = eval([stimtype '.dotwidth;']);
    wi = eval([stimtype '.stimwidth;']);
    hi = eval([stimtype '.stimheight;']);
    if ~isempty(dw)
        if sum(dw~=dw(1)) 
            eval([stimtype '.comment = [' stimtype '.comment ''NB different dotsizes were used'' ] ' ]);
        end
        if sum(wi~=wi(1))
            eval([stimtype '.comment = [' stimtype '.comment ''NB different stimulus widths were used'' ] ' ]);
        end
        if sum(hi~=hi(1)) 
            eval([stimtype '.comment = [' stimtype '.comment ''NB different stimulus heights were used'' ] ' ]);
        end
        
        eval([stimtype '.avgd.firingrate_mean = mean(' stimtype '.firingrate);']);
        eval([stimtype '.avgd.firingrate_SD = std(' stimtype '.firingrate);']);
        eval([stimtype '.avgd.nrep = length(' stimtype '.firingrate);']);
        eval([stimtype '.avgd.firingrate_sqmeansqrt = mean(sqrt(' stimtype '.firingrate))^2;']);
        eval([stimtype '.avgd.firingrate_meansqrt = mean(sqrt(' stimtype '.firingrate));']);
        eval([stimtype '.avgd.firingrate_SEM = ' stimtype '.avgd.firingrate_SD ./ sqrt(' stimtype '.avgd.nrep);']);
    end
end

% Get DTC
% Avergae over stimuli with identical disparities: record num rep, mean and SD
DTC.ndisp = length(DTC.disp_orthog);
if DTC.ndisp>0
    % If all stim params are equal except disparity:
    if sum(DTC.dotwidth~=DTC.dotwidth(1))~=0
        DTC.comment = [DTC.comment 'NB different dotsizes were used'];
    end
    if sum(DTC.disp_paral~=DTC.disp_paral(1))~=0
        DTC.comment = [DTC.comment 'NB different parallel disparities were used'];
    end
    if sum(DTC.stimwidth~=DTC.stimwidth(1))~=0
        DTC.comment = [DTC.comment 'NB different stimwidths were used'];
    end
    if sum(DTC.stimheight~=DTC.stimheight(1))~=0
        DTC.comment = [DTC.comment 'NB different stimheights were used'];
    end
    distinctdisp = DTC.disp_orthog(1);
    for j=2:DTC.ndisp
        if sum( distinctdisp == DTC.disp_orthog(j) )==0 % then we have not seen this disparity before
            distinctdisp = [distinctdisp DTC.disp_orthog(j)];
        end
    end
    
    DTC.avgd.disp_orthog = sort(distinctdisp);
    for j=1:length(distinctdisp)
        indx = find(DTC.disp_orthog == DTC.avgd.disp_orthog(j));
        DTC.avgd.firingrates{j} = DTC.firingrate(indx);
        DTC.avgd.firingrate_mean(j) = mean(DTC.firingrate(indx));
        DTC.avgd.firingrate_SD(j) = std(DTC.firingrate(indx));
        DTC.avgd.nrep(j) = length(indx);
        DTC.avgd.firingrate_sqmeansqrt(j) = mean(sqrt(DTC.firingrate(indx)))^2;
        DTC.avgd.firingrate_meansqrt(j) = mean(sqrt(DTC.firingrate(indx)));
    end
    DTC.avgd.firingrate_SEM = DTC.avgd.firingrate_SD ./ sqrt(DTC.avgd.nrep);        
    
    % Work out disparity discrimination index of Prince et al. 
    % NB include uncorrelated
    meansqrts = DTC.avgd.firingrate_meansqrt;
    nuncorr = length(uncorr.firingrate);
    nDTC = sum(DTC.avgd.nrep);
    m = nuncorr + nDTC;
    
    ressqrt = [];
    for j=1:length(DTC.avgd.firingrate_mean)
        ressqrt = [ressqrt (sqrt(DTC.avgd.firingrates{j}) - DTC.avgd.firingrate_meansqrt(j)).^2];
    end
    if nuncorr>0
        ressqrt = [ressqrt (uncorr.avgd.firingrate_meansqrt - sqrt(uncorr.firingrate)).^2 ];
        meansqrts = [meansqrts uncorr.avgd.firingrate_meansqrt];
    end
    Rmax = max(meansqrts);
    Rmin = min(meansqrts);
    n = length(meansqrts);
    RMSerror = sqrt(sum(ressqrt)/(m-n));
    RDSresponse.DDI = (Rmax - Rmin) / (Rmax - Rmin + 2*RMSerror) ;
    RDSresponse.SSQ = sum(ressqrt);
    RDSresponse.RMSerror = RMSerror;
    
end %if ndisp>0


RDSresponse.DTC.data = DTC;
RDSresponse.Lmonoc = Lmonoc;
RDSresponse.Rmonoc = Rmonoc;
RDSresponse.blank = blank;
RDSresponse.uncorr = uncorr;
RDSresponse.preperiod = preperiod;
