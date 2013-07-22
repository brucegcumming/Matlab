% Program to read in Bruce's rds.ODX.sum files.
function RDSresponse = ReadInRDS(RDSfile);
fid = fopen(RDSfile,'r');
if fid==-1
    % file could not be opened
    RDSresponse=0;
    return
end
RDSresponse.fileused = RDSfile;

while ~feof(fid)
    scrap = fgets(fid);
    % This file consists of lines with one word and one number
    [word,number] = strtok(upper(scrap));
    number = str2num(number);
    switch word
    case 'RDS_MONOC_LEFT'
        RDSresponse.Lmonoc = number;
    case 'RDS_MONOC_RIGHT'
        RDSresponse.Rmonoc = number;
    case 'RDS_BLANK_RESP'
        RDSresponse.blank = number;
    case 'RDS_UNCORR_RESP'
        RDSresponse.uncorr = number;
    case 'DISP_FREQ'
        RDSresponse.disparityfreq = number;
    end
end

fclose(fid)
