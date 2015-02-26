function DATA = CheckTrigTimes(DATA,varargin)
%DATA = CheckTrigTimes(DATA,varargin) Check DATA.trigtimes matches
%ClusterDetails
%

[DATA, Cd] = AllV.LoadClusterDetails(DATA);
Vall = AllV.mygetappdata(DATA,'Vall');

for j = 1:length(DATA.trigtimes)
    if ~isempty(DATA.trigtimes{j}) && isfield(Cd{j},'t')
        if length(DATA.trigtimes{j}) ~= length(Cd{j}.t)
             id = find(ismember(Vall.t,Cd{j}.t));
             DATA.trigtimes{j} = id;
        end        
    end
end