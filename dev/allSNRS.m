% find sessions with a spikes directory

clear all;close all;fclose all;
logfile='F:/Spike2/data/lem/SNRlog.txt';
SNRlog=fopen(logfile,'a');
if SNRlog==-1
    PrintMsg('could not open log file /n','%s');
    return
end
errcount=0;
sesscount=0;
for drive=1:2
    if drive==1
        topdir='F:/Spike2/data/lem/';
    else
        topdir='X:/smr/lem/';
    end
    cd(topdir);
    files=dir;
    for session=1:length(files)
        if strncmpi(files(session).name,'.',1)
            continue
        end
        if isdir(files(session).name)
            cd(files(session).name)
        else
            continue
        end
        sesscount=sesscount+1;
        if isdir('Spikes')
            PrintMsg(['Working on ',files(session).name ,'\n'],'%s');
            cd Spikes
            spkfiles=dir;
            for j=3:length(spkfiles)
                Es(j)=str2num(regexprep(spkfiles(j).name,['.*',files(session).name,'.p([0-9].*)t([0-9].*).mat'],'$2'));
            end
            for e=unique(Es(3:end))
                for p=1:24
                    try
                        load([files(session).name,'.p',num2str(p),'t',num2str(e),'.mat']);
                    catch
                        PrintMsg(SNRlog,['skipping ',files(session).name,' expt',' ',num2str(e),' because one or more spk files do not exist.' ,'\n'],'%s');
                        break
                    end
                    clusterfile=[topdir,files(session).name,'Expt',num2str(e),'ClusterTimes.mat'];
                    if exist(clusterfile)
                        auto=0;
                    else
                        auto=1;
                    end
                    try
                        res{sesscount}.stats{e,p}=getSNR(Spikes);
                    catch
                        errcount=errcount+1;
                        res{sesscount}.stats{e,p}={};
                        PrintMsg(SNRlog,['Problem with ',[files(session).name,'.p',num2str(p),'t',num2str(e),'.mat'],'. SNR could not be computed.' ,'\n'],'%s');
                        res{sesscount}.error{errcount}=lasterror;
                        rethrow(lasterror)
                        break
                    end
                    res{sesscount}.name=files(session).name;
                    res{sesscount}.drive=topdir(1:3);
                    res{sesscount}.auto{e,p}=auto;
                end
            end
            res{sesscount}.name=files(session).name;
            res{sesscount}.drive=topdir(1:3);
        end
        cd(topdir)
    end
end
   

%make matrix of SNR values
for sesscount=1:length(stats,1)
    for e=1:length(stats,2)
        for p=1:length(stats,3)
            if ~isempty(stats{sesscount,e,p}) & ~isempty(stats{sesscount,e,p}.SNR)
                SNRs(sesscount,e,p,:)=stats{sesscount,e,p}.SNR;
            else
                SNRs(sesscount,e,:,:)=NaN;
                break
            end
        end
    end
end

% find max SNR experiment for each session
SNRe=NaNmean(SNRs,3);
for sesscount=1:length(SNRe,1)
    for e=1:length(SNRe,2)
        if ~isnan(SNRe(sesscount,e+2))
            skipsum(sesscount,e)=sum(length(find(SNRe(sesscount,e,:)>3)),length(find(SNRe(sesscount,e+2,:)>3)));
        end
    end
end

%
for sesscount=1:length(skipsum,1)
    res{sesscount}.sessionSNRs=SNRs(sesscount,max(skipsum(sesscount,:)),:);
end



