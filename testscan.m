function Expts = testscan(name, varargin)
%Experiment with methods for reading psych files and other text data
%File format
% R%d  0 = wrong
%      1 = correct
%      2 = late/foul
%      3 = BadFix
%      >4 special lines with stimulus/expt info
%      11 Badfix caused by a microsaccade

% R%d xx=%f yy=%f time trialdur rwsize
% xx = value for expt 1
% yy = value for expt 2
%
% R = 0 = WRONG, 1 = CORRECT, 2 = FOUL Choice 3 = BAD-FIX
% R = 100 or 101 are 0,1 but in a correction loop
% 50,51 means saccade was not required
% R=109 indicates start of expt
% R9 Expt Start for Human Psych 
% R4 stimulus properties at start of expt, 
% R5 other stimulus properties
% R7 stimulus properties not in strict format - don't send these lines to
%                         textscan
% R27 Cancel Exp
% R10 = Expt End
% R8 = Expt Finished by verg


buildexpts = 1;
startday = 0;
DATA.verbose = 0;
DATA.plot.round = [0 0]; 
DATA.filename = name;
DATA.plot.mintrials = 10;
DATA.useallexpts = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'noexpts',8)
        buildexpts = 0;
    elseif strncmpi(varargin{j},'useallexpts',8)
        DATA.useallexpts = 1;
    end
    j = j+1;
end

tic;
txt = scanlines(name);
gid = find(~strncmp('R7',txt,2) & ~strncmp('testflag',txt,7));
xid = setdiff(1:length(txt),gid);

tic
for j = 1:length(gid)
    a = textscan(txt{gid(j)},'%s','delimiter',' ');
end
fprintf('textscan split %.4f\n',toc);

tic
a = textscan(char(txt(gid))','R%f %[^=]=%f %[^=]=%f %[^=]=%f %f %f %f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f');
fprintf('textscan took %.4f\n',toc)

tic
for j = 1:length(gid)
    a = textscan(txt{gid(j)},'R%f %[^=]=%f %[^=]=%f %[^=]=%f %f %f %f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f %[^=]=%f');
end
fprintf('textscan line by line took %.4f\n',toc);


tic
for j = 1:length(gid)
    a = split(txt{gid(j)});
end
fprintf('split took %.4f\n',toc);


tic
for j = 1:length(xid)
    s = split(txt{xid(j)});
end
fprintf('R7 split %d lines took %.4f\n',length(xid),toc);
