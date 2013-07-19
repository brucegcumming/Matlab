function cprintfig(varargin)
figure(2); plot(rand(1,100));
figure(1); plot(rand(1,100));
set(0,'currentfigure',2);
%figure(2); %This works. 
gcf
cprintf('blue','Figure now %d\n',gcf);
gcf
plot(rand(1,100),'r'); %Ths is in figure 1 !
