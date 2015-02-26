function DATA = BWPsychPlot(monkey,start,endpt,varargin)

% Plots a history of psychophysical performance given a range of dates,
% Looks in Z:/bgc/data/psych/ORBWRC for corresponding .mat files.
% start and endpt should be vectors of the form [a b c],
% where a,b,c are the "month", "day" and "year", as written in the filename.
% Note that this requires the inputs to match the strange date convention
% of the psych filenames, i.e. January = 0, and M.31.YYYY = (M+1).0.YYYY

if ispc
    drive='Z:/';
try cd(drive)
catch
    fprintf('check directory label. Problem finding network drive.')
    DATA=[];
    return
end
elseif ismac
    drive='/bgc/';
try 
    cd(drive)
catch
    drive='/Volumes/bgc/';
    try cd(drive)
    catch
    fprintf('check directory label. Problem finding network drive.')
    DATA=[];
    return
    end
end
end

resample=0;
rweight=0;
j=1;
latest=0;
while j <= length(varargin)
    if strncmpi(j,'rweight',2)
        rweight=1;
    elseif strncmpi(j,'resample',2)
        resample=1;
    elseif strncmpi(j,'latest')
        j=j+1;
        latest=varargin{j};
    end
    j=j+1;
end

if strcmpi(monkey,'lem')
    prefix='lem0000';
elseif strcmpi(monkey,'jbe')
    prefix='jbe000';
end
datelist = psychdatelist();
F=mat2str(datelist);
ststr=sprintf(';%s %s %s',num2str(start(1)),num2str(start(2)),num2str(start(3)));
enstr=sprintf(';%s %s %s',num2str(endpt(1)),num2str(endpt(2)),num2str(endpt(3)));
st=strfind(F,ststr)+1;
colon=strfind(F,';');
en=strfind(F,enstr)+1;
startit=find(colon<st);
endit=find(colon<en);
starti=startit(end)+1;
endi=endit(end)+1;
makelist=datelist(starti:endi,:);
[a b] = size(makelist);
figure;
k=0;
l=0;
for i=1:a
        days{i}=strcat(num2str(makelist(i,1)),'.',num2str(makelist(i,2)),'.',num2str(makelist(i,3)));
        files{i}=strcat(drive,'bgc/data/psych/ORBWRC/',monkey,'/',prefix,'.',days{i},'.mat');
    try 
        load(files{i});
        k=k+1;
        kernel{k}=Kernel;
        daystr{k}=strcat(num2str(makelist(i,1)),'.',num2str(makelist(i,2)),'.',num2str(makelist(i,3)));
        filename{k}=strcat(drive,'bgc/data/psych/ORBWRC/',monkey,'/',prefix,'.',daystr{k},'.mat');
    catch
        l=l+1;
        failstr{l}=strcat(num2str(makelist(i,1)),'.',num2str(makelist(i,2)),'.',num2str(makelist(i,3)));
        fprintf('\n no psychfile for %s',failstr{l});
        continue
    end
end
set(1,'Name',strcat(monkey,' kernels'));
if rweight
    BuildBWPsychRCa(kernel,'rweight');
else
   BuildBWPsychRCa(kernel); 
end

figure;
subplot(k,1,1);
set(gcf,'Name',strcat(monkey,' psychometric curves'));
for i=1:k
    subplot(k,1,i);
    ExptPsych(kernel{i},'shown');
end

figure;
subplot(a,1,1);
set(gcf,'Name',strcat(monkey,' Boot Distr. of Mn. Kernel Vector Angle'));
for i=1:k
    subplot(k,1,i);
    if isfield(kernel{i},'resampleors')
        F=find(kernel{i}.resampleors<0);
        kernel{i}.resampleors(F)=kernel{i}.resampleors(F)+180;
        hist(kernel{i}.resampleors,100);
        set(gca,'xlim',[0 180]);
        title(gca,daystr{i})
        if isfield(kernel{i},'recorddate')
            title(gca,[daystr{i} 'RECORDING DAY'])
        end
    else
        title(gca,strcat('No resamples made for ',daystr{i}));
    end
end

if resample
    for i=1:length(filename)
        a(i)=load(filename{i});
    end
    figure;
    if ~rweight
        Kernel=resampler(a);
    else
        Kernel=resampler(a,'rweight');
    end
        DATA.resample=Kernel;
end

DATA.filenames=filename;
DATA.kernels=kernel;


