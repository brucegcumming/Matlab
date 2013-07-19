function [DATA] = makenewORBWRCkernels(varargin)

%Searches through psych folders @
% 'Z:/bgc/data/psych/ORBWRC/' for any unconverted files.
%and makes into shortened '.mat' kernel files.  Prints addenda to logfile @
% 'Z:/bgc/data/psych/ORBWRC/kernellog.txt'
%must input desired monkey name(s) as argument inputs

if ispc
    drive='Z:/';
    log=fopen(strcat(drive,'bgc/data/psych/ORBWRC/kernellog.txt'),'a');
    if log==-1
        fprintf('check directory label. Problem finding network drive.')
    return
    end
elseif ismac
    drive='/bgc/';
    log=fopen(strcat(drive,'bgc/data/psych/ORBWRC/kernellog.txt'),'a');
    g=log;
    if g==-1
        drive='/Volumes/';
        log=fopen(strcat(drive,'bgc/bgc/data/psych/ORBWRC/kernellog.txt'),'a');
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

D=date;
fprintf(log,D);
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
    searchlist{j}=setdiff(searchilist{j},TreeFind(searchdir{j},'name','.*.mat'));
    L(j)=length(searchlist{j});
end
numstr=strfind(searchlist{1}{1},'RC/')+2;


%check which psych files are missing a corresponding .mat file (if no
%remake)
for j=1:2
    yes{j}=ones(1,length(searchlist{j}));
    counter=0;
    if ~remake
        outs{j}=cell(length(searchlist{j}));
        for i=1:L(j)
            A=strfind(searchlist{j}{i},'.mat');
            if isempty(A)
                continue
            else
                counter=counter+1;
            outs{j}{counter}=searchlist{j}{i}(1:A-1);
            end
        end
         for   i=1:counter
            if ismember(searchlist{j}{i},outs{j})
                yes{j}(i)=0;
            end
         end
    end
    need{j} = find(yes{j} == 1);
end

if isempty(need{1}) && isempty(need{2})
    PrintMsg(log,'There are no unmade kernels.');
    return
else
    lens=length(need{1})+length(need{2});
    PrintMsg(log,'%s kernels need to be made.',num2str(lens));
end   

%pass needed files to BuildBWPsychRCa and move subspace files
err=0;
for monkey=1:2        
    for i=1:length(need{monkey})
        if i>1
            PrintMsg(log,'Finished %s ',num2str(i-1));
            PrintMsg(log,' of %s files.',num2str(lens));
        end
        PrintMsg(log,'Working on %s . \n',searchlist{monkey}{need{monkey}(i)});
        try
            [a b] = BuildBWPsychRCa(searchlist{monkey}{need{monkey}(i)},'resample',1000,'save');
        catch
            err=err+1;
            PrintMsg(log,'%s could not be built.  Moving on.',searchlist{monkey}{need{monkey}(i)})
            DATA.err{monkey}{err}=searchlist{monkey}{need{monkey}(i)};
        end
        if isempty(a) && strncmpi(b.err,'Not ORBW',5)
            movedestination=strcat(savedir,monk{monkey},sscanf(searchlist{monkey}{need{monkey}(i)}(numstr+4:end),'%s'));
            movefile(searchlist{monkey}{need{monkey}(i)},movedestination);
            delete(strcat(searchlist{monkey}{need{monkey}(i)},'.mat'));
            PrintMsg(log,'Subspace file. Moved serial to %s and deleted .mat file.',movedestination);
        elseif isempty(a) && strncmpi(b.err,'Too Few',5)
            PrintMsg(log,'----> No kernel was made because too few trials. Keep .mat');
        elseif ~isempty(a)
            PrintMsg(log,' ----> Saved. \n');
        else
            PrintMsg(log,'unknown reason for empty .mat file. Kept in directory.');
        end
    end
end

return
    
    





