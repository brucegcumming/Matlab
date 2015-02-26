function Expt2Images(Expt, varargin)
%Rebuild RDS patterns used in an Experiment

seeds = [];
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'seeds',5)
        j = j+1;
        seeds = varargin{j};
    end
    j = j+1;
end

[a,filename] = fileparts(GetFileName(Expt));
imdir = ['/local/Images/' filename];
imdir = regexprep(imdir,'.[0-9]*$','');
options = unique({Expt.Trials.OptionCode});
optionflag = 'op=-cf-bw';
if sum(CellToMat(strfind(options,'+aa'))) == length(options)
    optionflag = [optionflag '+aa'];
end



txt{1} = 'mo=fore';
txt{end+1} = 'pausetimeout=0';
txt{end+1} = sprintf('st=rds\n\nfl=+pc+nc\n%s',optionflag);
txt{end+1} = sprintf('!savemovie=%s\n!renderoff2',imdir);
txt{end+1} = sprintf('!saveim12\nbc=0.5'); %defaults
S = Expt.Stimvals;
txt{end+1} = sprintf('uf=%s',filename);
f = {'px' 'py' 'vd' 'nf' 'wi' 'hi' 'dd' 'dw' 'bc' 'co' 'sl'};
for j = 1:length(f)
    if isfield(S,f{j})
        txt{end+1} = sprintf('%s=%.4f',f{j},S.(f{j}));
    end
end
if isfield(Expt.Stimvals,'mixac')
    txt{end+1} = sprintf('mixac=%.3f',Expt.Stimvals.mixac);
end

txt{end+1} = sprintf('uf=%s',filename);
if isempty(seeds)
    for j = 1:length(Expt.Trials)
        if isfield(Expt.Trials,'mixac')
            txt{end+1} = sprintf('mixac=%.3f',Expt.Trials(j).mixac);
        end
        if isfield(Expt.Trials,'dx')
            txt{end+1} = sprintf('dx=%.3f',Expt.Trials(j).dx);
        end
        for k = 1:length(f)
            if isfield(Expt.Trials,f{k})
                txt{end+1} = sprintf('%s=%.4f',f{k},Expt.Trials(j).(f{k}));
            end
        end
        se = Expt.Trials(j).se;
        id = Expt.Trials(j).id;
        txt{end+1} = sprintf('id=%.0f\nse=%.4f\n!onestim',id,se);
    end
else
    for k = 1:length(f)
       S.(f{k}) = GetEval(Expt,f{k});
       if ~isnan(S.(f{k}))
           txt{end+1} = sprintf('%s=%.4f',f{k},S.(f{k}));
       end
    end
    for j = 1:length(seeds)
        txt{end+1} = sprintf('se=%d\n!saveimage',seeds(j));
    end
end
    

fid = fopen('makerds.stm','w');
for j = 1:length(txt)
    fprintf(fid,'%s\n',txt{j});
end
fclose(fid);    


