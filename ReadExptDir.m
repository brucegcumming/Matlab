function AllExpts = ReadExptDir(name, varargin)
% Expts = ReadExptDir(name, varargin)
% Read all SPike2 .mat files in a directory and combines the Expts lists
% into one list
%
% Called by AplaySpkFile  if first argument is a directory

AllExpts = {};
%first sort numerically by suffix number
d = mydir([name '/*.mat']);
mnk =GetMonkeyName(name);
for j = 1:length(d)
    if ~isempty(regexp(d(j).name,[mnk '[M,\.,G,0-9]*.mat']))
        suffixes(j) = str2double(regexprep(d(j).filename,'.*[\.,A-z]([0-9]*).mat','$1'));        
    end
end
[a,b] =sort(suffixes);
sid = b(a> 0);


for j = 1:length(sid)
    fprintf('Reading %s\n',d(sid(j)).name);
    [Ex{j}, Expts] = APlaySpkFile(d(sid(j)).name,'nospikes','noerrs', varargin{:});
    AllExpts = {AllExpts{:} Expts{:}};
end

for j = 1:length(Ex)
    if isfield(Ex{j},'errs') && ~isempty(Ex{j}.errs)
        cprintf('blue','Errors for %s\n',d(sid(j)).name);
        for k = 1:length(Ex{j}.errs)
            cprintf('red','%s\n',Ex{j}.errs{k});
        end
    end
end

