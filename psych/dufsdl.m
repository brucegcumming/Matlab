if exist('diskprefix','var')
    prefix = [diskprefix '/bgc/data/dufus/psych/'];
else
    prefix = '\\lsr-bgc1/bgc/data/dufus/psych/';
end
alldata = [];
leftname = 'dufleft.mat';
rightname = 'dufright.mat';

prop = [];

for j =1:100;
  prop(j).xmin = -1;
  prop(j).xmax = 1;
  prop(j).nmin = 10;
  prop(j).sdk = 10;
  prop(j).sd = []; %force reading of only one sd
  prop(j).title = [];
  prop(j).forceread = 0;
end

nplot = 1;
alln = 1;
prop(alln).file = 'jul21';
prop(alln).skip = 1625;
prop(alln).last = 3500;
nplot = 1;

alln = alln+1;
prop(alln).file = 'jul21';
prop(alln).skip = 2;
prop(alln).last = 800;

alln = alln +1;
prop(alln).file = 'jul22';
prop(alln).skip = 200;
prop(alln).last = 1000;

alln = alln +1;
prop(alln).file = 'jul23';
prop(alln).skip = 100;
prop(alln).last = 1647;


prop(alln).title = 'No Delay';

alln = alln +1;
prop(alln).file = 'Aug06';
prop(alln).skip = 1;
prop(alln).last = 504;
prop(alln).sdk = 14;

alln = alln +1;
prop(alln).sdk = 14;
prop(alln).file = 'Aug07';
prop(alln).skip = 1;
prop(alln).last = 505;

alln = alln +1;
prop(alln).file = 'Aug08';
prop(alln).skip = 1;
prop(alln).last = 505;
prop(alln).sdk = 14;

alln = alln +1;
prop(alln).file = 'Aug11';
prop(alln).skip = 3;
prop(alln).last = 2013;
prop(alln).sdk = 14;

alln = alln +1;
prop(alln).file = 'Aug14';
prop(alln).skip = 1;
prop(alln).last = 1000;
prop(alln).sdk = 14;

alln = alln +1;
prop(alln).file = 'Aug15';
prop(alln).skip = 173;
prop(alln).last = 1291;
prop(alln).sdk = 14;

alln = alln +1;
prop(alln).file = 'Aug18';
prop(alln).sd = -1;
prop(alln).sdk = 14;

alln = alln +1;
prop(alln).file = 'Aug19';
prop(alln).skip = 778;
prop(alln).last = 1253;
prop(alln).sdk = 14;

alln = alln +1;
prop(alln).file = 'Aug20';
prop(alln).skip = 1;
prop(alln).last = 501;
prop(alln).sdk = 14;

alln = alln +1;
prop(alln).file = 'Aug20';
prop(alln).skip = 500;
prop(alln).last = 1117;
prop(alln).sdk = 14;

alln = alln +1;
prop(alln).file = 'Aug21';
prop(alln).skip = 16;
prop(alln).last = 539;
prop(alln).sdk = 14;

alln = alln +1;
prop(alln).file = 'Aug21';
prop(alln).skip = 539;
prop(alln).last = 1073;
prop(alln).sdk = 14;

alln = alln +1;
prop(alln).file = 'Aug22';
prop(alln).skip = 1;
prop(alln).last = 457;
prop(alln).sdk = 14;

alln = alln +1;
prop(alln).file = 'Aug25';
prop(alln).skip = 1;
prop(alln).last = 876;
prop(alln).sdk = 14;


alln = alln +1;
prop(alln).file = 'Aug21';
prop(alln).skip = 539;
prop(alln).last = 1073;
prop(alln).sdk = 14;


prop(alln).title  = 'Delay 14ms';

alln = alln +1;
prop(alln).file = 'Aug27';
prop(alln).skip = 432;
prop(alln).last = 1516;
prop(alln).sdk = 14;

%alln = alln +1;  % skip Aug27 - supicious
prop(alln).file = 'Sep02';
prop(alln).skip = 268;
prop(alln).last = 1391;
prop(alln).sdk = 14;

alln = alln +1;
prop(alln).file = 'Sep02';
prop(alln).skip = 1406;
prop(alln).last = 1999;
prop(alln).sdk = 14;

prop(alln).title  = 'Delay 28ms';


ranges = [];
ranges(1).min = -0.03;
ranges(1).max = 0.03;
ranges(1).sd = 14;
ranges(1).nmin = 100;

ranges(2).min = -0.03;
ranges(2).max = 0.03;
ranges(2).sd = -14;
ranges(2).nmin = 200;

ranges(3).min = -0.03;
ranges(3).max = 0.01;
ranges(3).sd = 0;
ranges(3).nmin = 200;

ranges(4).min = -0.07;
ranges(4).max = 0.07;
ranges(4).sd = -28;
ranges(4).nmin = 40;
side = 'L';

setnmin = 50;
makesummaries
summarize
suml = summ;
save(leftname,'suml');
xl = -0.03;
xu = 0.02;

