% Program to read in data on orientation tuning fom one of Bruce's "grating.OXM.fit" files
function ORTuning = ReadInORTC(ORgratingfile)

fid = fopen(ORgratingfile,'r');
if fid==-1
    % file could not be opened
    ORTuning=0;
    return
end
ORTuning.fit.fileused = ORgratingfile;
DEFINITIONS

while ~feof(fid)
    scrap = fgets(fid);
    
    for EYE=1:3
        if EYE==L
            tmp = strrep(scrap,'#FitL:','');
        elseif EYE==R
            tmp = strrep(scrap,'#FitR:','');
        elseif EYE==B
            tmp = strrep(scrap,'#FitB:','');
        end
        
        if ~strcmp(tmp,scrap) % then the line began with "#Fit"
            % Read data both for fit with baseline, and for fit without baseline
            for j=1:2
                if j==1
                    CONSTRAINT = NOBASE;
                elseif j==2
                    % Need to read in the next line
                    scrap = fgets(fid);
                    CONSTRAINT = WITHBASE;
                    tmp = strrep(scrap,'#Base:','');
                end
                fitparams = sscanf(tmp,'%f');
                if length(fitparams)~=4
                    disp('I was expecting 4 numbres here ... something wrong')
                    stop
                end
                if ~isempty(findstr(tmp,'Lin')) % then this is the linear fit
                    FIT=LIN;
                elseif ~isempty(findstr(tmp,'Log')) % then this is a log fit
                    FIT=LOG;
                end
                % Convert so that thetapref really is the preferred OR, not its log
                if FIT==LIN
                    ORTuning.fit.params(EYE,FIT,CONSTRAINT).thetapref = fitparams(1);
                elseif FIT==LOG
                    ORTuning.fit.params(EYE,FIT,CONSTRAINT).thetapref = exp(fitparams(1));
                end
                ORTuning.fit.params(EYE,FIT,CONSTRAINT).SD = fitparams(2);
                ORTuning.fit.params(EYE,FIT,CONSTRAINT).peak = fitparams(3);
                ORTuning.fit.params(EYE,FIT,CONSTRAINT).base = fitparams(4);
                % Remove everything up to the 'R' that indicates Residual:
                indx=findstr(tmp,'R')+1;
                tmp2 = tmp(indx:end);
                % Read off residual and percentage of variance
                tmp2 = sscanf(tmp2,'%f %f');
                ORTuning.fit.params(EYE,FIT,CONSTRAINT).RSS = tmp2(1);
                ORTuning.fit.params(EYE,FIT,CONSTRAINT).varexplnd = tmp2(2);
                if CONSTRAINT == WITHBASE
                    % Remove everything up to the '%' that indicates % var explained:
                    indx=findstr(tmp,'%')+1;
                    tmp2 = tmp(indx:end);
                    tmp2 = sscanf(tmp2,'%f');
                    % check to see whether baseline was worth including
                    ORTuning.fit.params(EYE,FIT,CONSTRAINT).pbaseline = tmp2;
                end

            end % end loop over base/no base
        end % end if-loop checking that first line begins with "#Fit:
    end % end loop over eyes

    tmp = strrep(scrap,'Preferred Eye','');
    if ~strcmp(tmp,scrap) % then this is the line recording the preferred eye & which fit was best
        if strcmp(tmp(end-1),'R')
            ORTuning.fit.DomEye = R;
        elseif strcmp(tmp(end-1),'L')
            ORTuning.fit.DomEye = L;
        elseif strcmp(tmp(end-1),'B')
            ORTuning.fit.DomEye = B;      
        else
            disp('I was expecting this to be labelled either L or R or B')
            stop
        end
    end
    
    tmp = strrep(scrap,'FitAreaL','');
    if ~strcmp(tmp,scrap) % then this is the line recording area under fitted curve
        ORTuning.fit.areabestfitORTC(L) = sscanf(tmp,'%e');
    end
    
    tmp = strrep(scrap,'FitAreaR','');
    if ~strcmp(tmp,scrap) % then this is the line recording area under fitted curve
        ORTuning.fit.areabestfitORTC(R) = sscanf(tmp,'%e');
    end
    
    tmp = strrep(scrap,'DataArea','');
    if ~strcmp(tmp,scrap) % then this is the line recording area under data
        if ~isempty(findstr(tmp,'Left'))
            ORTuning.data.areaORTC(L) = sscanf(tmp,'%e');
        elseif ~isempty(findstr(tmp,'Right'))
            ORTuning.data.areaORTC(R) = sscanf(tmp,'%e');
        elseif ~isempty(findstr(tmp,'Binoc'))
            ORTuning.data.areaORTC(B) = sscanf(tmp,'%e');
        else
            disp('I was expecting this to be labelled either L or R or B')
            stop
        end 
    end
        
end % end while loopgoing through file
    
    
% For each eye, work out whether lin or log fit was best, and whether it needed a baseline
[neyes,nfits,nconst] = size(ORTuning.fit.params);
ORTuning.fit.BestFitType = [0 0 0];
ORTuning.fit.BestConstraint = [0 0 0];
for EYE=1:neyes
    if ~isempty(ORTuning.fit.params(EYE,FIT,WITHBASE).varexplnd)
        for FIT=[LIN] % only linear fits
            varexplained(FIT) = max([ORTuning.fit.params(EYE,FIT,WITHBASE).varexplnd ORTuning.fit.params(EYE,FIT,NOBASE).varexplnd]);
        end
        [mx,bestfittype] = max(varexplained);
        ORTuning.fit.BestFitType(EYE) = bestfittype;
        % Record whether baseline was worth including as a free param
        if ORTuning.fit.params(EYE,bestfittype,WITHBASE).pbaseline < 0.05
            ORTuning.fit.BestConstraint(EYE) = WITHBASE;
        else
            ORTuning.fit.BestConstraint(EYE) = NOBASE;
        end
    end
end
fclose(fid)