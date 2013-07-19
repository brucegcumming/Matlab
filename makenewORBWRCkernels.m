function [DATA] = makenewORBWRCkernels(varargin)

%Searches through psych folders @
% 'Z:/bgc/data/psych/ORBWRC/' for any unconverted files.
%and makes into shortened '.mat' kernel files.  Prints addenda to logfile @
% 'Z:/bgc/data/psych/ORBWRC/Kernel Logs/kernellog.txt'
%must input desired monkey name(s) as argument inputs

D=date;
t=now;
if ispc
    drive='Z:/';
    log=fopen(strcat(drive,'bgc/data/psych/ORBWRC/Kernel Logs/kernellog-',D,num2str(t),'.txt'),'a');
    if log==-1
        fprintf('check directory label. Problem finding network drive.')
    return
    end
elseif ismac
    drive='/bgc/';
    log=fopen(strcat(drive,'bgc/data/psych/ORBWRC/Kernel Logs/kernellog-',D,num2str(t),'.txt'),'a');
    g=log;
    if g==-1
        drive='/Volumes/';
        log=fopen(strcat(drive,'bgc/data/psych/ORBWRC/Kernel Logs/kernellog-',D,num2str(t),'.txt'),'a');
        if log==-1;
            g=-1;
        else
            g=1;
        end
        if g==-1;
            fprintf('check directory label. Problem finding network drive.')
            return
        end
    end
end
path(strcat(drive,'bgc/matlab'),path);
path(strcat(drive,'bgc/matlab/dev'),path);
g=1;
remake=0;
monk={};
prefix={};

while g<=length(varargin)
    if strcmpi(varargin{g},'remake')
        remake=1;
    elseif strcmpi(varargin{g},'lem')
        monk{g}='lem';
        prefix=[prefix,'lem0000'];
    elseif strcmpi(varargin{g},'jbe')
        monk{g}='jbe';
        prefix=[prefix,'jbe000'];
    end
    g=g+1;
end
if isempty(monk)
    PrintMsg('No monkeys selected!')
    return
elseif isempty(monk{1})
    monk=monk{2:end};
    prefix=prefix{2:end};
end

%set directories of interest
savedir = strcat(drive,'bgc/data/psych/ORSubspace/');
for i=1:2
searchdir{i} = strcat(drive,'bgc/data/psych/ORBWRC/',monk{i});
end

%find psych files
for j=1:2
    searchylist{j}=TreeFind(searchdir{j},'name',prefix{j});
    k=0;
    for i=1:length(searchylist{j})
        if strfind(searchylist{j}{i},'#')
        else
            k=k+1;
            searchilist{j}{k}=searchylist{j}{i};
        end
    end
    outlist{j}=TreeFind(searchdir{j},'name','.*.mat');
    searchlist{j}=setdiff(searchilist{j},outlist{j});
end
numstr=strfind(searchlist{1}{1},'RC/')+2;
DATA=struct('result',cell(1,2),'filenames',searchlist);

%check which psych files are missing a corresponding .mat file (if no
%remake)
for j=1:2
    yes{j}=zeros(1,length(searchlist{j}));
    if ~remake
        for i=1:length(outlist{j})
            A=strfind(outlist{j}{i},'.mat');
            outs{j}{i}=outlist{j}{i}(1:A-1);
        end
    end
    yes{j} = ismember(searchlist{j},outs{j});
    need{j} = find(yes{j} == 0);
end

if isempty(need{1}) && isempty(need{2})
    PrintMsg(log,'There are no unmade kernels.');
    return
else
    lens=length(need{1})+length(need{2});
    PrintMsg(log,'%s kernels need to be made.',num2str(lens));
end   

%pass needed files to BuildBWPsychRCa and move subspace files
for monkey=1:2        
    DATA(monkey).result=cell(1,length(need{monkey}));
    for i=1:length(need{monkey})
        tim=tic;
        save(strcat(drive,'bgc/data/psych/ORBWRC/Kernel Logs/DATAlog',D,num2str(t),'.mat'),'DATA');
        fprintf(log,'%s',t);
        if i>1
            PrintMsg(log,'Finished %s of %s files. \n',num2str(i-1),num2str(lens));
        end
        PrintMsg(log,'Working on %s . \n',searchlist{monkey}{need{monkey}(i)});
        try
            [a b] = BuildBWPsychRCa(searchlist{monkey}{need{monkey}(i)},'resample',1000,'save');
        catch
            error=lasterror;
            PrintMsg(log,error.message);
            result=PrintMsg(log,'%s could not be built.  Moving on.\n',searchlist{monkey}{need{monkey}(i)})
            DATA(monkey).result{i}=error;
            save(strcat(searchlist{monkey}{need{monkey}(i)},'.mat'),'error');
            tom=toc(tim);
            PrintMsg(log,'Took %s seconds.',num2str(round(tom)));
            continue
        end
        if isempty(a) && strncmpi(b.err,'Not ORBW',5)
            movedestination=strcat(savedir,monk{monkey},sscanf(searchlist{monkey}{need{monkey}(i)}(numstr+4:end),'%s'));
            movefile(searchlist{monkey}{need{monkey}(i)},movedestination);
            delete(strcat(searchlist{monkey}{need{monkey}(i)},'.mat'));
            result=PrintMsg(log,'Subspace file. Moved serial to %s and deleted .mat file.\n',movedestination);
        elseif isempty(a) && strncmpi(b.err,'Too Few',5)
            result=PrintMsg(log,'----> No kernel was made because too few trials. Keep .mat');
        elseif ~isempty(a)
            result=PrintMsg(log,' ----> Saved. \n');
        else
            result = PrintMsg(log,'unknown reason for empty .mat file. Kept in directory.');
        end
        DATA(monkey).result{i}=result;
        tom=toc(tim);
        PrintMsg(log,'Took %s seconds.',num2str(round(tom)));
    end
end
return




