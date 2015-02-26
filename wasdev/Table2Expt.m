function matexp = Table2Expt(name, varargin)


expname = [];
basedir = [];
lines = [];    
nrpt = 0;
matexp.binocstrs = {};

j = 1;
while j <= length(varargin)
    if strcmp(varargin{j},'basedir',7)
        j = j+1;
        basedir = varargin{j};
    elseif strcmp(varargin{j},'expname',7)
        j = j+1;
        expname = varargin{j};
    end
    j = j+1;
end

if ischar(name)
    lines = scanlines(name);
    diroot = fileparts(name);
end

started = 0;
ns = 0;
teststim = [];
if ~isempty(lines)
    for j = 1:length(lines)
        if isempty(deblank(lines{j}))
        elseif started
            f = split(lines{j});
            if length(f) > length(expvars)
                fprintf('line %d too many  fields\n',j);
            elseif length(f) == length(expvars)
                ns = ns+1;
                if f{1}(1) == '*'
                    teststim(end+1) = ns;
                    f{1} = f{1}(2:end);
                end
                for k = 1:length(f)
                    if f{k}(1) == '$' %special syntax for predefined values
                        nv = sscanf(f{k}(2:end),'%d');
                        if nv && nv <= length(values)
                            f{k} = values{nv};
                        end
                    end
                    AllS(ns).(expvars{k}) = f{k};
                    if strcmp(expvars{k},'st')
                        exvals(ns,k) = StimulusName(f{k});
                    elseif strcmp(expvars{k},'op')
                        exvals(ns,k) = 0;
                    else
                    exvals(ns,k) = str2num(f{k});
                    end
                end
            elseif length(f) > 1  %error.
            end
        else
            id = strfind(lines{j},'=');
            if ~isempty(id)
                val = lines{j}(id(1)+1:end);
            else
                val = '';
            end
            if lines{j}(1) == ':'
                expvars = split(lines{j}(2:end));
                started = 1;
            elseif strncmp(lines{j},'values:',7)
                f = split(lines{j});
                values = f(2:end);                
            elseif strncmp(lines{j},'basedir',7)
                if isempty(basedir) && ~isempty(val)
                    basedir = val;
                end
            elseif strncmp(lines{j},'expname',7)
                if isempty(expname) && ~isempty(val)
                    expname = val;
                end
            else
                matexp.binocstrs{end+1} = lines{j};
            end        
        end
    end
end

ns = length(AllS);

if isempty(basedir) && ~isempty(diroot)
    basedir=diroot;
    fprintf('Using %s as dir\n',basedir);
end

for j = 1:ns
    s{j} = WriteStim(basedir, j-1, AllS(j), exvals(j,:));
    fprintf('%d: %s\n',(j-1),s{j});
end

if nrpt == 0  %set nr automatically to get ~ 80 trials
nrpt = round(80./ns);
end

if isempty(teststim)
    stimorder = repmat([1:ns]-1,1,nrpt);
    stimorder = stimorder(randperm(length(stimorder)));
else
    stimorder = repmat(teststim-1,1,nrpt);
end
f = fields(AllS);

fid = fopen([basedir '/stimorder'],'w');
fprintf(fid,'expvars %s',f{1});
for j = 2:length(f)
    fprintf(fid,',%s',f{j});
end
fprintf(fid,'\n');
if ~isempty(expname)
    fprintf(fid,'expname=%s\n',expname);    
    matexp.expname = expname;
else
fprintf(fid,'expname=NotSet\n');    
end
fprintf(fid,'%s\n',sprintf('%d ',stimorder));
fclose(fid);
fprintf('%d stim * %d repeats\n',ns,nrpt);
matexp.AllS = AllS;
matexp.stimdir=basedir;
matexp.exvals = exvals;

function s = WriteStim(basedir, stimno, S, exvals)

stimname = sprintf('%s/stim%d',basedir,stimno);
fid = fopen(stimname,'w');
f = fields(S);
exstr = [];
s = [];
for j = 1:length(f)
    x = S.(f{j});
    str = sprintf('%s=%s',f{j},x);
    fprintf(fid,'%s\n',str);
    s = [s str ' '];
end
exstr = sprintf(' %.2f',exvals);
fprintf(fid,'manexpvals%d%s\n',stimno,exstr);

fclose(fid);

