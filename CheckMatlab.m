function CheckMatlab(name)
%CheckMatlab(name) Count number of lines in functions for a matlab source file 
%Analyses size of source (.m) files
if ~exist(name,'file')
    name = which(name);
end
fid = fopen(name);
if fid
    a = textscan(fid,'%s','delimiter','\n');
    txt = a{1};
    fcid = regexp(txt,'\s*function');
    for j = length(fcid):-1:1
        fc(j) = sum(fcid{j});
    end
    lines = find(fc);
    lens = diff(lines);
    lens(length(lines)) = length(txt)-lines(end);
    [a,b] = sort(lens);
    for j = 1:length(lens)
        fprintf('%d %s\n',lens(b(j)),txt{lines(b(j))});
    end
    hist(lens);
    title(sprintf('%d functions, mean %.1f lines',length(lens),mean(lens)));
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