function timerfna(tim, varargin)
DATA = get(findobj('Tag',get(tim,'Tag')),'UserData');
%DATA = get(findobj('Tag','Combiner'),'UserData');
if DATA.state.autolist
tic;
cmb.combine('relist');
%toc
end

