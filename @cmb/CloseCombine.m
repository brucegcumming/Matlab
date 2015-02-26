function CloseCombine(a,b, varargin)

if nargin == 0  %manual call
end
DATA = GetDataFromFig(a);

j = 1;
args = {};
while j <= length(varargin)
    if strncmpi(varargin{j},'combinespecial',11)
    else
        args = {args{:} varargin{j}};
    end
    j = j+1;
end

if DATA.state.online && isfield(DATA,'spkcache') %save classifications so far
    outname = [DATA.name '/SpikeClassification.mat'];
    if ~isempty(DATA.spkcache)
        spkdata = DATA.spkcache;
%not sure this is worth the time - have to reload spk files anyway...        
%        save(outname,'spkdata');
    end
end
cmb.DoLayout(DATA,'savelast');
cmb.DoConfig(DATA,'savelast');
if isfield(DATA,'timerobj')
    stop(DATA.timerobj);
end
if isfield(DATA,'tag')
    CloseChildren(DATA.toplevel, args{:});
    names = fieldnames(DATA.tag);
    for j = 1:length(names)
        CloseTag(DATA.tag.(names{j}));
    end
end

