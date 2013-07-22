% Program to read in Bruce's rds.ODX.fit files.
function DTC = Read_rds_ODX_fit(DTCfile);
fid = fopen(DTCfile,'r');
if fid==-1
    % file could not be opened
    DTC=0;
    return
end
DTC.fit.fileused = DTCfile;

while ~feof(fid)
    scrap = fgets(fid);
    % This file consists of lines with one word and one number
    [word,number] = strtok(upper(scrap));
    if ~isempty(word)
        switch word
        case 'CORRELATED'        
            [word,number] = strtok(upper(number));
            switch word
            case 'FIT:'
                params = str2num(number);
                DTC.fit.baseline = params(2);
                DTC.fit.amp = params(3);
                DTC.fit.SF = params(4);
                DTC.fit.phase = params(5);
                DTC.fit.SD = params(6);
                DTC.fit.offset = params(7);        
            end
        end
    end
end

fclose(fid)