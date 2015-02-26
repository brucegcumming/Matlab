ENDSTIM=3;
FRAMESIGNAL=5;
BADFIX=11;
STARTTRIAL=12;
STARTEXPT=1;
ENDEXPT=2;
STARTSTIM=6;
WURTZLATE=7;
WURTZOK=8;
WURTZPREM=9;
WURTZOKW=20;
ENDTRIAL=16;
CANCELEXPT=19;
BWISREADY=4;

STOREBIT=16;
PSYCHBIT = 2^15;

stimnames = {'none' ,	'gabor',	'rds' ,	'grating',	'bar',	'circle',...
	'rectangle','test',	'square',	  'probe',	  '2grating',  'cylinder',...
	  'corrug',	'sqcorrug',	'twobar',	'rls', 'annulus', 'rdssine', 'nsines', 'rlssine',...
	  'radial', 'image', 'checker'};
  
  LMONOC=1;
  RMONOC=2;

 IUNCORR = -1005;
 IBLANK = -1009;
 ILEFTMONOC = -1001;
 IRIGHTMONOC = -1002;
 ISIGNALFRAME = -1010;
 IBROADBAND = -1011;
 
 j = 1;
 CodeNames.Codes{j} = 'or'; CodeNames.Label{j} = 'Ori'; j = j+1;
 CodeNames.Codes{j} = 'jv'; CodeNames.Label{j} = 'Speed'; j = j+1;
 CodeNames.Codes{j} = 'sz'; CodeNames.Label{j} = 'Size'; j = j+1;
 CodeNames.Codes{j} = 'sf'; CodeNames.Label{j} = 'SF'; j = j+1;
 CodeNames.Codes{j} = 'tf'; CodeNames.Label{j} = 'TF'; j = j+1;
 CodeNames.Codes{j} = 'me'; CodeNames.Label{j} = 'Eye'; j = j+1;
 CodeNames.Codes{j} = 'dw'; CodeNames.Label{j} = 'DotSz'; j = j+1;
 CodeNames.Codes{j} = 'ed'; CodeNames.Label{j} = 'ed'; j = j+1;
 CodeNames.Codes{j} = 'Dw'; CodeNames.Label{j} = 'Dw'; j = j+1;
  
 