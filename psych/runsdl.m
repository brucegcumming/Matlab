if exist('diskprefix','var')
    prefix = [diskprefix '/bgc/data/rufus/psych/'];
else
    prefix = '\\lsr-bgc1/bgc/data/rufus/psych/';
end

leftname = 'rufleft.mat';
rightname = 'rufright.mat';

alldata = [];
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
prop(alln).file = 'Jul04';
prop(alln).xmin = -0.05;
prop(alln).xmax = 0.05;
prop(alln).sdk = 14;
alln = alln-1;

alln = alln+1;
prop(alln).file = 'Jul08';
prop(alln).xmin = -0.05;
prop(alln).xmax = 0.05;
prop(alln).sdk = 14;
alln = alln-1;

alln = alln+1;
prop(alln).file = 'Jul09';
prop(alln).xmin = -0.05;
prop(alln).xmax = 0.05;
prop(alln).sdk = 14;
alln = alln-1;

alln = alln+1;
prop(alln).file = 'Nov17';
prop(alln).xmin = -0.05;
prop(alln).xmax = 0.05;
prop(alln).sdk = 14;

prop(alln).title = 'No Delay';

alln = alln+1;
prop(alln).file = 'Aug19';
prop(alln).last = 978;
prop(alln).sdk = 14;

alln = alln+1;
prop(alln).file = 'Aug20';
prop(alln).last = 1424;
prop(alln).sdk = 14;

alln = alln+1;
prop(alln).file = 'Aug21';
prop(alln).skip = 24;
prop(alln).last = 943;
prop(alln).sdk = 14;

alln = alln+1;
prop(alln).file = 'Aug22';
prop(alln).skip = 571;
prop(alln).last = 1344;
prop(alln).sdk = 14;

alln = alln+1;
prop(alln).file = 'Aug22';
prop(alln).skip = 1937;
prop(alln).last = 2243;
prop(alln).title = '28 ms';
prop(alln).sdk = 14;

alln = alln+1;
prop(alln).sd = -1;
prop(alln).file = 'Jul22';
prop(alln).last = 1752;
prop(alln).xmin = -0.1;
prop(alln).sdk = 14;
prop(alln).xmax = 0.06;

alln = alln+1;
prop(alln).file = 'Jul18_1';
prop(alln).last = 1752;
prop(alln).xmin = -0.1;
prop(alln).xmax = 0.1;
prop(alln).sd = -1;
prop(alln).sdk = 14;

alln = alln+1;
prop(alln).file = 'Jul29';
prop(alln).last = 1752;
prop(alln).xmin = -0.1;
prop(alln).xmax = 0.1;
prop(alln).sdk = 14;

alln = alln+1;
prop(alln).file = 'Jul30_1';
prop(alln).last = 143;
prop(alln).xmin = -0.1;
prop(alln).xmax = 0.1;
prop(alln).sdk = 14;

prop(alln).title = 'SD 14 ms';

alln = alln+1;
prop(alln).file = 'Aug06';
prop(alln).skip = 100;
prop(alln).xmin = -0.05;
prop(alln).xmax = 0.05;

alln = alln+1;
prop(alln).file = 'Nov19';
prop(alln).skip = 0;
prop(alln).xmin = -0.05;
prop(alln).xmax = 0.05;

alln = alln+1;
prop(alln).file = 'Aug07';
prop(alln).skip = 20;
prop(alln).last = 807;
prop(alln).xmin = -0.05;
prop(alln).xmax = 0.05;

alln = alln+1;
prop(alln).file = 'Aug08_10';
prop(alln).last = 601;
prop(alln).xmin = -0.05;
prop(alln).xmax = 0.05;

alln = alln+1;
prop(alln).file = 'Aug11';
prop(alln).skip = 1215;
prop(alln).xmin = -0.05;
prop(alln).xmax = 0.05;

prop(alln).title = 'SD 10ms';

alln = alln+1;
prop(alln).file = 'Aug12';
prop(alln).last = 999;
prop(alln).xmin = -0.05;
prop(alln).xmax = 0.05;

alln = alln+1;
prop(alln).file = 'Aug13';
prop(alln).last = 2005;
prop(alln).xmin = -0.05;
prop(alln).xmax = 0.05;

%Aug 13 breaks code. So skip by not incrementing alln
%alln = alln+1;
prop(alln).file = 'Aug14';
prop(alln).last = 1913;
prop(alln).xmin = -0.05;
prop(alln).xmax = 0.05;

alln = alln+1;
prop(alln).file = 'Aug15';
prop(alln).skip = 1102;
prop(alln).last = 1604;
prop(alln).xmin = -0.05;
prop(alln).xmax = 0.05;

alln = alln+1;
prop(alln).file = 'Aug18';
prop(alln).last = 1789;
prop(alln).xmin = -0.05;
prop(alln).xmax = 0.05;




makesummaries;

ranges = [];
ranges(1).sd = 14;
ranges(1).min = -0.02;
ranges(1).max = 0.02;
ranges(1).nmin = 70;

side = 'L';
ylim = [0.005 0.2];
summarize;
suml = summ;
save(leftname,'suml');
