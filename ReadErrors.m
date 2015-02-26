function ReadErrors(name, varargin)
%print out errors for Data Folders
savefile = 0;

j = 1;
while j <= length(varargin)
    if strcmp(varargin{j},'save')
        savefile = 1;
    end
    j = j+1;
end


if iscellstr(name)
    for j = 1:length(name)
        ReadErrors(name{j}, varargin{:});
    end
    return;
end

if isdir(name)
    errfile = [name '/Errors.mat'];
    if exist(errfile)
        load(errfile);
        missing = 0;
    else
        missing = 1;
        Errors = {};
    end
    d = dir([name '/backup/Errors*.mat']);
    for j = 1:length(d)
        ename = [name '/backup/' d(j).name];
        X = load(ename);
        Errors = {Errors{:} X.Errors{:}};
    end
    
    ShowErrors(Errors);
    if (savefile || missing) && ~isempty(Errors)
        fprintf('Saving %s\n',errfile);
        save(errfile,'Errors');
    end
end


function ShowErrors(E)
for j = 1:length(E)
    fprintf('%s\n',E{j}.s);
end