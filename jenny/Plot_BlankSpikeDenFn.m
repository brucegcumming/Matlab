function Plot_BlankSpikeDenFn(results)
% Plot spike density function during blank screen stimulus
spkfile = strrep(results.RDSresponse.DTC.data.fileused,'.trl','.spk');
% Assume cluster is 1...
cluster = 1;
% ... unless told otherwise
indx = findstr(spkfile,'_');
if ~isempty(indx)
    cluster = sscanf(spkfile((indx+1):end),'%d');
end

[tspike,stimonsetindx,stimduration] = ReadBlankSpikeTimes(spkfile,cluster);
nstimruns = length(stimonsetindx);

% Work out spike density functions 
mintimetoplot = 0.050;
maxtimetoplot = max(stimduration)+0.050;

time = [ mintimetoplot : (maxtimetoplot-mintimetoplot)/1e4 : maxtimetoplot ];
spikedenfun = zeros(size(time));
tSD = 1e-3; % Width of Gaussian convolved with each spike

for jstim=1:nstimruns
    spktimes = tspike{jstim};
    nspks(jstim) = length(spktimes);
    for jspk=1:nspks(jstim)
        spikedenfun = spikedenfun + exp(-(time-spktimes(jspk)).^2./2/tSD^2);
    end % next spike
end % next stimulus
 
% normalise so that area under curve divided by total time gives mean firing rate
% Need to correct for how many runs went into it
spikedenfun = spikedenfun / sqrt(2*pi) / tSD  / nstimruns;

save BlankSpikeDenFn

figure('numbertitle','off','name',results.filebase,'pos',[78         737        1002         373])
% mark on screen rate
FT = trapz(time,spikedenfun .* exp(2*pi*i*72.*time));
phase = angle(FT);
framedur = 1/72;
t=mintimetoplot-framedur-framedur/4+framedur*phase/2*pi;
mx = 1.1*max(spikedenfun);
pale = [1 1 1]*0.7;
while t<maxtimetoplot
    patch([t t+framedur/2 t+framedur/2 t],[0 0 mx mx],pale,'edgecol',pale)
    t = t + framedur;
end
hold on
area(time,spikedenfun,'facecolor',[1 1 1]*0.,'edgecolor','k')
xlabel('time (seconds)','fontsize',14)
ylabel('instantaneous firing rate','fontsize',14)
axis tight
set(gca,'xlim',[mintimetoplot maxtimetoplot])
title(strrep(results.filebase,'.0',''),'fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tspike_myclst,stimonsetindx,stimduration] = ReadBlankSpikeTimes(file,cluster)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In each file a line of the form:
% Stim 24571121, 2.005 dx=0.360 or=170 sf=4.000 tf=8.000 true=8.168
% tells you the disparity (dx=...) of the stimulus, and the start time.
% Lines of the form
% S 4211 0.15 -44 4357 4289 0.270 0.002
% describe minsaccades that happpendduring the stimulus.
% and lines like:
% E es 24591167
% describe events that occured during the stimulus, all of which you
% can ignore. The Event es indicates the end of the stimulus, but this
% is always a fixed period so you don't really need to worry.
% All lines beginning with a digit indicat the time of a spike. The fist digit
% is the clster classification of the spike, the second number is the spike 
% time. clster 0 is multiunit junk and should be ignored. In most cases
% what you want is clster 1, but in a few clster 2 is the disparity  tuned
% cell. THe list below says which clster you want in each file.
% In this routine, I extract only those spikes which occurred where the stimulus was a blank screen

fid = fopen(file,'r');
% Go through and look for a line beginning with "Type"
scrap='';
while isempty(strmatch('Type',scrap))
    scrap = fgetl(fid);
end

% Go through and look for a line beginning with "Stim"
stimrun=0;
while ~feof(fid)
   scrap='';
   while isempty(strmatch('Stim',scrap)) & ~feof(fid)
       scrap = fgetl(fid);
   end
   if ~feof(fid) & ... now this is where I select only those Stim lines specifying a blank
           isempty(findstr(scrap,'lm')) & isempty(findstr(scrap,'rm')) & isempty(findstr(scrap,'uc')) & isempty(findstr(scrap,'se')) & ~isempty(findstr(scrap,'cn'))
       stimrun=stimrun+1;       
       scrap = strrep(strrep(scrap,'Stim',''),',','')
       tmp = sscanf(scrap,'%d %f');
       stimonsetindx(stimrun) = tmp(1);
       stimduration(stimrun) = tmp(2);
       
       tspike = ReadData(fid,cluster);
       % Convert to seconds from 1e4s:
       tspike = tspike.*1e-4
       % Only record spikes that occur after 50ms hvae elapsed since stim onset, and before 50ms after stim offset
       tspike = tspike(tspike>0.050 & tspike <  stimduration(stimrun)+0.050);
       
       % Record these
       tspike_myclst{stimrun} = tspike
   end
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  tspike = ReadData(fid,cluster)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scrap='';
ndata=0;
tspike = [];
% Continue until next line marked End is found
while isempty(strmatch('End',scrap)) & ~feof(fid)
   %position = ftell(fid);
   scrap = fgets(fid);
   thisclster = str2num(sscanf(scrap,'%s',1));
   % Look for lines beginning with the desired cluster index:
   if ~isempty(thisclster)
       if thisclster==cluster
           ndata=ndata+1;
           % I want to read the second element of the line
           tmp = sscanf(scrap,'%e',[1 2]);
           tspike(ndata) = tmp(2);
       end
   end
end

% Reset the file positon indicator so I get the Stim line again hen I next read it
%if ~feof(fid); fseek(fid,position,'bof'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [dx,or,sf,tf,true] = ExtractParams(filine)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Have reached a line beginning with "Stim". Proceed to extract disparity ('dx') etc.
% i) Extract disparity: dx=???
tmp = filine;
scrap='';
j=0;while isempty(strmatch('dx=',scrap)) & j<1e3; j=j+1;   scrap = sscanf(tmp,'%s',1); tmp = strrep(tmp,scrap,''); end
dx = str2num(strrep(scrap,'dx=',''));
% ii) Extract orientation: or=???
tmp = filine;
j=0;while isempty(strmatch('or=',scrap)) & j<1e3; j=j+1;   scrap = sscanf(tmp,'%s',1); tmp = strrep(tmp,scrap,''); end
or = str2num(strrep(scrap,'or=',''));
% iii) Extract spatial frequency: sf=???
tmp = filine;
j=0;while isempty(strmatch('sf=',scrap)) & j<1e3; j=j+1;   scrap = sscanf(tmp,'%s',1); tmp = strrep(tmp,scrap,''); end
sf = str2num(strrep(scrap,'sf=',''));
% iv) Extract timefreq: tf=???
tmp = filine;
j=0;while isempty(strmatch('tf=',scrap)) & j<1e3; j=j+1;   scrap = sscanf(tmp,'%s',1); tmp = strrep(tmp,scrap,''); end
tf = str2num(strrep(scrap,'tf=',''));
% iv) Extract ??: true=???
tmp = filine;
j=0;while isempty(strmatch('true=',scrap)) & j<1e3; j=j+1;   scrap = sscanf(tmp,'%s',1); tmp = strrep(tmp,scrap,''); end
true = str2num(strrep(scrap,'true=',''));

