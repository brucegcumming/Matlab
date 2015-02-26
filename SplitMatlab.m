function SplitMatlab(name, outname)
%SplitMatlab(name, @class) Splits a .m file into separate files for each function
%and puts them into a class folder. The name for this must include the "@"
if ~exist(name,'file')
    name = which(name);
end
fid = fopen(name);
if fid
    a = textscan(fid,'%s','delimiter','\n');
    txt = a{1};
    fcid = regexp(txt,'^\s*function');
    for j = length(fcid):-1:1
        fc(j) = sum(fcid{j});
    end
    lines = find(fc);
    prefix = outname(2:end);
    outfile = [outname '/' prefix '.m'];
    fid = fopen(outfile,'w');
    fprintf(fid,'classdef %s\nproperties\n Version = ''1.1''\nend\n',prefix);
    fprintf(fid, 'methods (Static)\n');
    for j = 1:length(lines)
        flin = txt{lines(j)};
        fprintf('%s\n',flin);
        fprintf(fid,'%s\n',regexprep(flin,'function\s+',''));
        flin = regexprep(flin,'function .* =\s*','');
        flin = regexprep(flin,'function\s+','');
        fname = regexprep(flin,'(.*','');
        fnames{j} = fname;
    end
    fprintf(fid,'\nend\nend\n');
    fclose(fid);
    for j = 1:length(txt)
        if ~ismember(j,lines)
            for k = 1:length(fnames)
                if strfind(txt{j},fnames{k})
                    if isempty(strfind(txt{j},[fnames{k} '('])) && isempty(strfind(txt{j},['@' fnames{k}]))
                        cprintf('blue','%s:%sNotreplaced\n',txt{j},fnames{k});
                    end
                end
                txt{j} = strrep(txt{j},[fnames{k} '('],[ prefix '.' fnames{k} '(']);
                txt{j} = strrep(txt{j},['@' fnames{k}],['@' prefix '.' fnames{k}]);
            end
        end
    end
    for j = 1:length(lines)
        outfile = [outname '/' fnames{j} '.m'];
        fprintf('Writing %s\n',outfile);
        fid = fopen(outfile,'w');
        if j < length(lines)
            endline = lines(j+1)-1;
        else
            endline = length(txt);
        end
        for k = lines(j):endline
            fprintf(fid,'%s\n',txt{k});
        end
        fclose(fid);
    end
    v = ver('MATLAB');
    version = sscanf(v.Version(3:end),'%d');
    if version  > 7
        a = checkcode(name);
        for j = 1:length(a)
            if ~isempty(strfind(a(j).message,'function')) && ~isempty(strfind(a(j).message,'unused'))
                fprintf('%s\n',a(j).message);
            end
        end
    end
end

function b = noused(a)

b = sqrt(a);