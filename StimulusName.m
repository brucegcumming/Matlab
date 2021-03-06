function name = StimulusName(st)
%Convert stimulus name string to num or vice versa
stimnames = {'none' ,	'gabor',	'rds' ,	'grating',	'bar',	'circle',...
	'rectangle','test',	'square',	  'probe',	  '2grating',  'cylinder',...
	  'corrug',	'sqcorrug',	'twobar',	'rls', 'annulus', 'rdssine', 'nsines', 'rlssine',...
	  'radial', 'image', 'checker'};
  
if ischar(st)
    name = find(strcmp(st, stimnames));
    if ~isempty(name)
        name = name-1;
    else
        name = NaN;
    end
else
    name = stimnames{st+1};
end
