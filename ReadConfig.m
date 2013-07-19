function DATA = ReadConfig(DATA, configfile, varargin)
%Read a file that sets the configuation of GUI parameters
%See also SaveConfig

%     if ~exist(configfile)
%         return;
%     end

    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'print',4)
            fprintf('Getting settings from %s\n',configfile);
        end
        j = j+1;
    end
    

    try
        fid = fopen(configfile,'r');
    catch
         fprintf('Cant read %s\n',configfile);
         return;
    end
    if fid > 0
    a = textscan(fid,'%s','delimiter','\n');
    for j = 1:length(a{1})
        try
            eval(a{1}{j});
        catch
            fprintf('Bad Syntax %s\n',a{1}{j});
        end
    end
    fclose(fid);
    end