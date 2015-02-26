function LFP = FixLFPFile(name, varargin)
%
%LFP = FixLFPFile(name, varargin)
% given name of .mat file for expt (NOT the LFP file)
% loads relevant LFP file, and Fixes for Mains noise, if necessary.

logfid = 0;
LFP = [];
mainsname = 'Ch24';
mainsch = 24;

j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'logfid',6)
        j = j+1;
        logfid = varargin{j};
    end
    j = j+1;
end

if exist(name,'file');
load(name);
else
    return;
end

lfpfile = strrep(name,'.mat','A.lfp.mat');
fixfile = strrep(name,'.mat','.lfp.mat');
% with 8 channels, everything is in one file - there is no 'A.lfp.mat'
if ~exist(lfpfile,'file') && exist(fixfile,'file')
    lfpfile = fixfile;
    forcefix = 1;
end

    vars = who;
    for j = 1:length(vars)
        if ~isempty(regexp(vars{j},'Ch[0-9][0-9]*'))
            eval(['ch = ' vars{j} ';']);
            if strncmpi(ch.title,'Mains',5)
                mainsname = vars{j};
                mainsch = ch;
            end
        end
    end

if exist(lfpfile,'file') && (~exist(fixfile,'file') || forcefix) ...
        && exist(mainsname,'var')
    load(lfpfile);
    if isfield(LFP.Header,'MainsNoise')
        fprintf('%s Already fixed\n',lfpfile);
        if logfid
            fprintf(logfid,'%s Already fixed\n',lfpfile);
        end
    else
        LFP = FixLFPMains(LFP,mainsch.times .* 10000);
        LFP.Header.amps = LFPGains(LFP);
        a = LFP.Header.amps ./ max(LFP.Header.amps);
        nch = sum(LFP.Header.chanlist > 0);
        if std(a(find(a > 0.1))) > 0.2 && nch <= 8  %% 8 channel probe tends to have mixed LFP gain
            LFP.Header.needscale = 1;
        else
            LFP.Header.needscale = 0;
        end
        save(fixfile,'LFP');
    end
end

