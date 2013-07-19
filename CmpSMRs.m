function [matchfiles, details] = CmpSMRs(monkey,varargin)
%
% CmpSMRs(monkey) finds files for monkey on C:\Spike6\data and on Z:\smr
% and lists differences
%
% CmpSMRs(monkey,'localdir', lpath) searches for local files in lpath,
% instead of C:\Spike6\data
% CmpSMRs(monkey,'netdir', npath) searches for local files in npath
%
% CmpSMRs(monkey,'suffix', suffix) searches for files with extensinon suffix instead of .smr

localdir = ['C:\Spike6\data\' monkey];
netdir = ['Z:\smr\' monkey];
matchfiles = {};
suffix = '.smr$';
suffixes{1} = suffix;
details = [];

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'localdir',8)
        j = j+1;
        localdir = varargin{j};
    elseif strncmpi(varargin{j},'netdir',8)
        j = j+1;
        netdir = varargin{j};
    elseif strncmpi(varargin{j},'suffix',8)
        j = j+1;
        suffix = varargin{j};
    end
    j = j+1;
end
[locals, lsizes, ldates] = TreeFind(localdir,'name',suffix);
[nets, nsizes, ndates] = TreeFind(netdir,'name',suffix);
for j = 1:length(nets)
    [ndirs{j}, netnames{j}] = fileparts(nets{j});
end

for j = 1:length(locals)
    [ldir, name] = fileparts(locals{j});
    id = strmatch(name,netnames);
    if length(id) > 0
        matches(j) = id(1);
        if nsizes(id(1)) == lsizes(j) && ndates(id(1)) >= ldates(j)
            result(j) = 0;
        elseif nsizes(id(1)) == lsizes(j) && ndates(id(1)) < ldates(j)
            result(j) = 0;
        else
            result(j) = 1;
        end
    else
        result(j) = 2;
    end
end


n = 0;
id = find(result == 0);
for j = id
   fprintf('%s matches %s: %.1fMb\n',locals{j},nets{matches(j)},lsizes(j)./1048576);
   n = n+1;
   matchfiles{n} = locals{j};
end
fid = id;
id = find(result == 1);
n = 0;
for j = id
            fprintf('Size Mismatch %s is %.1fMb, %s is %.1fMb\n',locals{j},lsizes(j)./1048576,nets{matches(j)},nsizes(matches(j))./1048576);
   n = n+1;
   details.sizewrong{n} = locals{j};
end
mid = id;
n = 0;
id = find(result == 2);
for j = id
        fprintf('No match for %s (%.1fMb)\n',locals{j},lsizes(j)./1048576);
   n = n+1;
   details.unmatch{n} = locals{j};
end

for j = 1:n
    [a,b] = fileparts(details.unmatch{j});
    allb{j} = regexprep(b,'\..*','');
    allb{j} = regexprep(allb{j},'A$','');
end
if n > 0
details.unmatchfile = unique(allb);
end
n = 0;
id = find(result == 3);
for j = id
        fprintf('%s newer than %s\n',locals{j},nets{matches(j)});
   n = n+1;
   details.newer{n} = locals{j};
end

fprintf('Totals (Gb) Match %.2f, Mismatch %.2f, Unmatch %.2f\n',sum(lsizes(fid))./(2^30),sum(lsizes(mid))./(2^30),sum(lsizes(id))./2^30);