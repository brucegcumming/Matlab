% find sessions with a spikes directory

clear all;close all;fclose all;
logfile='F:/Spike2/data/lem/manSNRlog.txt';
SNRlog=fopen(logfile,'a');
if SNRlog==-1
    PrintMsg('could not open log file /n','%s');
    return
end
errcount=0;
sesscount=0;
cd 'X:/smr/lem/';
sessions={'M195' 'M203' 'M211'};
for session=1:length(sessions)
    cd(sessions{session})
    sesscount=sesscount+1;
    PrintMsg(['Working on ',sessions{session} ,'\n'],'%s');
    load CellList
    cd Spikes
    cells=unique(CellList(:));
    a=find(cells>0);
    cells=cells(a);
    for cell=cells'
        ind = find(CellList==cell);
        [E P cluster] = ind2sub(size(CellList),ind);
        l=length(E);
        for count=1:l
            for p=P'
                e=E(count);
                try
                    load([sessions{session},'.p',num2str(p),'t',num2str(e),'.mat']);
                catch
                    res{sesscount}.stats{cell}{count}={};
                    PrintMsg(SNRlog,['skipping ',sessions{session},' expt',' ',num2str(e),' because spk file ''' ,[sessions{session},'.p',num2str(p),'t',num2str(e),'.mat'] ,''' does not exist.' ,'\n'],'%s');
                    break
                end
                try
                    res{sesscount}.stats{cell}{count}=getSNR(Spikes,'cluster',cluster(count));
                catch
                    errcount=errcount+1;
                    res{sesscount}.stats{cell}{count}={};
                    PrintMsg(SNRlog,['Problem with ',[sessions{session},'.p',num2str(p),'t',num2str(e),'.mat'],'. SNR could not be computed.' ,'\n'],'%s');
                    res{sesscount}.error{errcount}=lasterror;
                    res{sesscount}.error{errcount}.session=sessions{session};
                    res{sesscount}.error{errcount}.ep = [e p];
                    res{sesscount}.error{errcount}.cellno = cell;
                    break
                end
            end
        end
        res{sesscount}.name=sessions{session};
    end
cd 'X:/smr/lem/';
end

%    
% %make matrix of SNR values
% for sesscount=1:length(stats,1)
%     for e=1:length(stats,2)
%         for p=1:length(stats,3)
%             if ~isempty(stats{sesscount,e,p}) & ~isempty(stats{sesscount,e,p}.SNR)
%                 SNRs(sesscount,e,p,:)=stats{sesscount,e,p}.SNR;
%             else
%                 SNRs(sesscount,e,:,:)=NaN;
%                 break
%             end
%         end
%     end
% end
% 
% % find max SNR experiment for each session
% SNRe=NaNmean(SNRs,3);
% for sesscount=1:length(SNRe,1)
%     for e=1:length(SNRe,2)
%         if ~isnan(SNRe(sesscount,e+2))
%             skipsum(sesscount,e)=sum(length(find(SNRe(sesscount,e,:)>3)),length(find(SNRe(sesscount,e+2,:)>3)));
%         end
%     end
% end
% 
% %
% for sesscount=1:length(skipsum,1)
%     res{sesscount}.sessionSNRs=SNRs(sesscount,max(skipsum(sesscount,:)),:);
% end



