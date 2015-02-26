function MakeProbeMenu(DATA, sm, varargin)

if nargin == 1 || isempty(sm)
    sm = findobj(DATA.toplevel,'type','uimenu','tag','ProbeSwitchMenu');
    if isempty(sm)
        return;
    end
end
delete(get(sm,'children'));
pchars = ['1':'9' '0' 'a':'z'];
Clusters = AllV.mygetappdata(DATA, 'Clusters');
for j = 1:DATA.allnprobes
    if j < 10
        uimenu(sm,'Label',['&' num2str(j)] ,'Callback',{@AllV.ChangeProbe, j});
    elseif j <= length(pchars)
        uimenu(sm,'Label',[num2str(j) ' (&'  pchars(j) ')'] ,'Callback',{@AllV.ChangeProbe, j});
    else
        uimenu(sm,'Label', num2str(j) ,'Callback',{@AllV.ChangeProbe, j});
    end
    if length(Clusters) >= j && isfield(Clusters{j},'trigset')
        p = j + 0.1;
        uimenu(sm,'Label',[num2str(j) ' Set' num2str(2)] ,'Callback',{@AllV.ChangeProbe, p});
    end
end
