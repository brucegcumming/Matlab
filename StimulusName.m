function name = StimulusName(st)

stimnames = {'none' ,	'gabor',	'rds' ,	'grating',	'bar',	'circle',...
	'rectangle','test',	'square',	  'probe',	  '2grating',  'cylinder',...
	  'corrug',	'sqcorrug',	'twobar',	'rls', 'annulus', 'rdssine', 'nsines', 'rlssine',...
	  'radial', 'image', 'checker'};
  
if ischar(st)
    name = strcmp(st, stimnames);
    if ~isempty(name)
        name = name-1;
    else
        name = NaN;
    end
else
    name = stimnames(st+1);
end
