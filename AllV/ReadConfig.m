function DATA = ReadConfig(DATA, configfile, varargin)
%Read a file that sets the configuation of GUI parameters
%See also SaveConfig

%     if ~exist(configfile)
%         return;
%     end
usegui = 0;
chooseifnone = 1; %if file doesn't exist, use gui to choose

t = 'Read Configuration';

    j = 1;
    while j <= length(varargin)
        if strncmpi(varargin{j},'print',4)
            fprintf('Getting settings from %s\n',configfile);
        elseif strncmpi(varargin{j},'choose',4)
            usegui = 1;
        elseif strncmpi(varargin{j},'nochoose',4)
            chooseifnone = 0;
        end
        j = j+1;
    end
    
 if ~exist(configfile) && chooseifnone
     t = sprintf('File does not exist %s',configfile);
     usegui = 1;
 end
 
    if usegui
        [outname, path] = uigetfile(configfile,t);
        if ischar(outname)
            configfile = [path outname];
        else
            return;
        end
    end

    try
        fid = fopen(configfile,'r');
    catch
         cprintf('errors','Cant read %s\n',configfile);
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
 