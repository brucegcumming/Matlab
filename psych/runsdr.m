
if exist('diskprefix','var')
    prefix = [diskprefix '/bgc/data/rufus/psych/'];
else
    prefix = '\\lsr-bgc1/bgc/data/rufus/psych/';
end
endalldata = [];
prop = [];

leftname = 'rufleft.mat';
rightname = 'rufright.mat';

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
prop(alln).file = 'Aug27';
prop(alln).skip = 1;
prop(alln).last = 418;

alln = alln+1;
prop(alln).file = 'Aug25';
prop(alln).skip = 64;
prop(alln).last = 662;
prop(alln).title = 'No Delay';

alln = alln+1;
prop(alln).sd = [];
prop(alln).file = 'Aug27';
prop(alln).skip = 428;
prop(alln).last = 1120;
prop(alln).xmin = -0.025;
prop(alln).xmax = 0.014;

alln = alln+1;
prop(alln).file = 'Aug28';
prop(alln).skip = 179;
prop(alln).last = 888;

alln = alln+1;
prop(alln).file = 'Aug28';
prop(alln).skip = 1432;
prop(alln).last = 1825;
prop(alln).sd = -1;
prop(alln).title = 'SD 10 ms';

alln = alln+1;
prop(alln).file = 'Aug29';
prop(alln).skip = 1295;
prop(alln).last = 1825;
prop(alln).sdk = 14;

alln = alln+1;
prop(alln).file = 'Aug30';
prop(alln).skip = 5;
prop(alln).last = 603;
prop(alln).sdk = 14;

alln = alln+1;
prop(alln).file = 'Aug30';
prop(alln).skip = 603;
prop(alln).last = 1107;
prop(alln).sdk = 14;
prop(alln).title = 'SD 14ms';

alln = alln+1;
prop(alln).file = 'Aug31';
prop(alln).skip = 15;
prop(alln).last = 517;
%prop(alln).xmin = -0.05;
%prop(alln).xmax = 0.05;

alln = alln+1;
prop(alln).file = 'Aug31';
prop(alln).skip = 1006;
prop(alln).last = 1205;

alln = alln+1;
prop(alln).file = 'Sep01';
prop(alln).skip = 84;
prop(alln).last = 494;
prop(alln).sdk = 14;

alln = alln+1;
prop(alln).file = 'Sep01';
prop(alln).skip = 896;
prop(alln).last = 1297;
prop(alln).title = 'SD 20ms';
prop(alln).sdk = 14;

alln = alln+1;
prop(alln).file = 'Sep01';
prop(alln).skip = 1344;
prop(alln).last = 1835;
prop(alln).sdk = 14;

alln = alln+1;
prop(alln).file = 'Sep02';
prop(alln).skip = 1263;
prop(alln).last = 0;
prop(alln).xmin = -0.035;
prop(alln).xmax = 0.0205;
prop(alln).sdk = 14;

alln = alln+1;
prop(alln).file = 'Sep02';
prop(alln).skip = 448;
prop(alln).last = 1054;
prop(alln).sdk = 14;
prop(alln).title = 'SD 28ms';


makesummaries;

ranges(1).sd = 10;
ranges(1).min = -0.025;
ranges(1).max = 0.014;
ranges(1).nmin = 10;
ranges(2).sd = -20;
ranges(2).min = -0.1;
ranges(2).max = 0.1;
ranges(2).nmin = 10;
ranges(3).sd = 0;
ranges(3).min = -0.02;
ranges(3).max = 0.02;
ranges(3).nmin = 20;
ranges(4).sd = 28;
ranges(4).min = -0.16;
ranges(4).max = 0.1;
ranges(4).nmin = 50;
ranges(5).sd = 28;
ranges(5).min = -0.1;
ranges(5).max = 0.1;
ranges(5).nmin = 10;

side = 'R';
summarize;
sumr = summ;
save('rightall.mat','sumr');
