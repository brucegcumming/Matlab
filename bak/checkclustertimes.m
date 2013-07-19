function [suspicious] = checktimeranges(directory)
tic;
%checks for cluster files that have suspiciously short numbers of trials.
%plots a histogram with time ranges for all manually cut clusters.
%the data in the histogram is also saved as Times, a 1xn cell array, where
%n is the number of experiments.  Times{n} is a 1x24 vector, with a time
%range for each probe.  NaN means it was autocut and therefore the field
%restricttimerange does not exist.

Times=cell(1,500);

for i=1:500
    Times{i}=cell(1,24);
end

for i=1:500
    for j=1:24
        Times{i}{j}=NaN;
    end
end



cd(directory)

files=dir;

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
              
             gg= 'caught1';
       
          end
      end
          
      
      else
          
       gg= 'caught2';
          
    end
    
    catch
        
   gg=' caught3';
   
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

figure;hist(histinput);

[H G]=find(histinput>0);F=histinput(G);[gg hh]=sort(F);

[t tt]=find(gg<5);
t=max(tt);

list(1,:)=sort(G(hh(1:t)));

list(2,:)=gg(1:t);

for i=1:length(list(1,:))
    suspicious(i).probe=rem(list(1,i),24);
    suspicious(i).expt=files((list(1,i)-suspicious(i).probe)/24+1).name;
    suspicious(i).duration=list(2,i);
end

toc
