function [suspicious] = checktimeranges(directory)
%given a directory, such as 'F:/Spike2/data/lem/M195'
%checktimeranges checks for cluster files that have suspiciously short numbers of trials.
%output: 'suspicious' is a structured array with fields giving the experiment and probe with suspiciously short time
%ranges, and shows the time range.
tic;
cd(directory)
files=dir;
Times=cell(1,length(files));
for i=1:length(files)
    Times{i}=cell(1,24);
end
for i=1:length(files)
    for j=1:24
        Times{i}{j}=NaN;
    end
end
for i=1:length(files)
      try
        if files(i).name(6) == 'C'
            clear Clusters
            clear ClusterDetails
            load(files(i).name);  
        elseif files(i).name(7) == 'C' 
            clear Clusters
            clear ClusterDetails
            load(files(i).name);
            for j=1:24
               try
                   Times{i}{j}=Clusters{j}.restricttimerange(2)-Clusters{j}.restricttimerange(1);
               catch
               end
            end
        else        
        end
      catch
      end
end
k=0;
histinput=zeros(1,length(files)*24);
for i=1:length(files)
    for j=1:24
        k=k+1;
        histinput(k)=Times{i}{j};
    end
end
%1. plots a histogram with time ranges for all manually cut clusters.
figure;hist(histinput);
set(gcf,'Name','distribution of time ranges for manually restricted experiments');
[H G]=find(histinput>0);F=histinput(G);[gg hh]=sort(F);
[t tt]=find(gg<5);
if isempty(t);
suspicious='There are no suspiciously short experiments';
else
t=max(tt);
list(1,:)=sort(G(hh(1:t)));
list(2,:)=gg(1:t);
for i=1:length(list(1,:))
    suspicious(i).probe=rem(list(1,i),24);
    suspicious(i).expt=files((list(1,i)-suspicious(i).probe)/24+1).name;
    suspicious(i).duration=list(2,i);
end
toc
end
