prop = [];

if exist('diskprefix','var')
    prefix = [diskprefix '/bgc/data/dufus/psych/'];
else
    prefix = '\\lsr-bgc1/bgc/data/dufus/psych/';
end
alldata = [];
side = 'R';
leftname = 'dufleft.mat';
rightname = 'dufright.mat';

for j =1:100;
  prop(j).xmin = -1;
  prop(j).xmax = 1;
  prop(j).nmin = 10;
  prop(j).sdk = 10;
  prop(j).sd = []; %force reading of only one sd
  prop(j).title = [];
  prop(j).skip = [];
  prop(j).last = [];
  prop(j).forceread = 0;
end

nplot = 1;
alln = 1;
prop(alln).file = 'Oct15';
prop(alln).xmin = -0.09;
prop(alln).xmax = 0.09;

%alln = alln+1;  % discard Oct 15 - was bad
prop(alln).file = 'Oct16';
prop(alln).xmin = -0.09;
prop(alln).xmax = 0.09;

alln = alln +1;
prop(alln).file = 'Oct17';
prop(alln).xmin = -0.09;
prop(alln).xmax = 0.09;
prop(alln).last = 1009;
alln = alln - 1;

alln = alln +1;
prop(alln).file = 'Nov06';
prop(alln).xmin = -0.09;
prop(alln).xmax = 0.09;
prop(alln).last = 1020;
alln = alln - 1;

alln = alln +1;
alln = 1;
prop(alln).file = 'Nov20';
prop(alln).skip = 100;
prop(alln).last = 600;
prop(alln).xmax = 0.08;
prop(alln).xmin = -0.08;

alln = alln +1;
prop(alln).file = 'Nov21';
prop(alln).skip = 100;
prop(alln).last = 0;
prop(alln).xmax = 0.04;
prop(alln).xmin = -0.04;


prop(alln).title = 'No Delay';
alln = alln +1;
prop(alln).file = 'Oct21';
prop(alln).xmin = -0.09;
prop(alln).xmax = 0.09;
prop(alln).xmax = 0.025;

alln = alln +1;
prop(alln).file = 'Oct22';
prop(alln).xmin = -0.09;
prop(alln).xmax = 0.025;
prop(alln).skip = 1;
prop(alln).last = 1006;

alln = alln +1;
prop(alln).file = 'Oct23';
prop(alln).xmin = -0.09;
prop(alln).xmax = 0.025;
prop(alln).skip = 1;

alln = alln +1;
prop(alln).file = 'Oct24';
prop(alln).xmin = -0.03;
prop(alln).xmax = 0.03;
prop(alln).skip = 1;
prop(alln).last = 500;

alln = alln +1;
prop(alln).file = 'Oct28';
prop(alln).xmin = -0.09;
prop(alln).xmax = 0.025;
prop(alln).title = '10ms';

alln = alln +1;
prop(alln).file = 'Oct31';
prop(alln).xmin = -0.05
prop(alln).xmax = 0.05;
prop(alln).sdk = 14;

alln = alln +1;
prop(alln).file = 'Nov04';
prop(alln).xmin = -0.05
prop(alln).xmax = 0.05;
prop(alln).sdk = 14;
prop(alln).skip = 166;

alln = alln +1;
prop(alln).file = 'Nov05';
prop(alln).xmin = -0.05
prop(alln).xmax = 0.05;
prop(alln).sdk = 14;

alln = alln +1;
prop(alln).file = 'Dec10';
prop(alln).xmin = -0.05
prop(alln).xmax = 0.05;
prop(alln).sdk = 14;
prop(alln).last = 1000;

alln = alln +1;
prop(alln).file = 'Dec11';
prop(alln).xmin = -0.05
prop(alln).xmax = 0.05;
prop(alln).sdk = 14;
prop(alln).last = 1000;

alln = alln +1;
prop(alln).file = 'Dec12';
prop(alln).xmin = -0.05
prop(alln).xmax = 0.05;
prop(alln).sdk = 14;
prop(alln).last = 0;

alln = alln +1;
prop(alln).file = 'Dec17';
prop(alln).xmin = -0.05
prop(alln).xmax = 0.05;
prop(alln).sdk = 14;
prop(alln).last = 2000;

prop(alln).title = '14 ms';

alln = alln +1;
prop(alln).file = 'Dec03';
prop(alln).xmin = -0.05
prop(alln).xmax = 0.05;
prop(alln).sdk = 10;

alln = alln +1;
prop(alln).file = 'Dec04';
prop(alln).xmin = -0.05
prop(alln).xmax = 0.02;
prop(alln).sdk = 10;

alln = alln +1;
prop(alln).file = 'Dec05';
prop(alln).xmin = -0.05
prop(alln).xmax = 0.05;
prop(alln).sdk = 10;

alln = alln +1;
prop(alln).file = 'Dec09';
prop(alln).xmin = -0.05
prop(alln).xmax = 0.05;
prop(alln).sdk = 10;
prop(alln).title = '20 ms';

alln = alln +1;
prop(alln).file = 'Dec18';
prop(alln).xmin = -0.15
prop(alln).xmax = 0.15;
prop(alln).sdk = 14;

alln = alln +1;
prop(alln).file = 'Dec19';
prop(alln).xmin = -0.05
prop(alln).xmax = 0.05;
prop(alln).sdk = 14;
alln = alln -1; % skip Dec 19 - Bad

alln = alln +1;
prop(alln).file = 'Dec23';
prop(alln).xmin = -0.05
prop(alln).xmax = 0.05;
prop(alln).sdk = 14;

alln = alln +1;
prop(alln).file = 'Dec24';
prop(alln).xmin = -0.05
prop(alln).xmax = 0.05;
prop(alln).sdk = 14;
prop(alln).title = '28 ms';




ranges =[];
ranges(1).sd = -10;
ranges(1).min = -0.03;
ranges(1).max = 0.01;
ranges(1).nmin = 50;

ranges(2).sd = 10;
ranges(2).min = -0.03;
ranges(2).max = 0.01;
ranges(2).nmin = 50;

ranges(3).sd = 0;
ranges(3).min = -0.03;
ranges(3).max = 0.01;
ranges(3).nmin = 50;

ranges(4).sd = -20;
ranges(4).min = -0.05;
ranges(4).max = 0.05;
ranges(4).nmin = 50;

ranges(5).sd = -14;
ranges(5).min = -0.05;
ranges(5).max = 0.05;
ranges(5).nmin = 100;

ranges(6).sd = 14;
ranges(6).min = -0.03;
ranges(6).max = 0.03;
ranges(6).nmin = 200;

ranges(7).sd = 20;
ranges(7).min = -0.05;
ranges(7).max = 0.05;
ranges(7).nmin = 50;

ranges(8).sd = 28;
ranges(8).min = -0.15;
ranges(8).max = 0.15;
ranges(8).nmin = 50;

ranges(9).sd = -28;
ranges(9).min = -0.08;
ranges(9).max = 0.08;
ranges(9).nmin = 50;

makesummaries;
summarize;
sumr = summ;
save(rightname,'sumr');


xl = -0.12;
xu = 0.12;
