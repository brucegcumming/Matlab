 function  [true, each] =  NeedTemplateForCluster(C, all) needothermeans = 0; if ~isfield(C,'space') %%No cluster     true = 0;     each = [];     return; end each(1) = AllV.WhichPlotType(C, 1); if sum(ismember(C.space,[17 18]))     needothermeans = 1; end if all && isfield(C,'next') for j = 1:length(C.next)         each(j+1) = AllV.WhichPlotType(C,j+1);         if isfield(C.next{j},'space') && sum(ismember(C.next{j}.space,[17 18]))             needothermeans = 1;         end end end true = sum(ismember(each,[3 4])) > 0; if needothermeans      true = 2; end  %Find flags/states set in the GUI andcopy them to a new DATA struct%