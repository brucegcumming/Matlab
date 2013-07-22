function DATA = MergExcludedTrials(DATA)

if isfield(DATA,'CellDetails') && isfield(DATA.CellDetails,'excludetrials')
    
    for j = 1:length(DATA.CellDetails.exptids)
        eid = floor(DATA.CellDetails.exptids(j));
        for k = 1:size(DATA.CellDetails.excludetrials,2)
            for c = 1:size(DATA.CellDetails.excludetrials,3)
                if ~isempty(DATA.CellDetails.excludetrials{j,k,c})
                    DATA.AllClusters{eid}(k).excludetrialids{c} = union(DATA.AllClusters{eid}(k).excludetrialids{c},DATA.CellDetails.excludetrials{j,k,c});
                end
            end
        end
    end
end