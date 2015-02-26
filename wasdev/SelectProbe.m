function selected = SelectProbe(probelist, selected, varargin)
%selected = SelectProbe(probelist, selected)
%popup with check boxes for selecting probes
%selected = SelectProbe(probelist, selected,  'single')
%       returns as soon as one checkbox is hit

fontsize = 14;
singleselection = 0;
Array = [];
if isfield(probelist,'X')
    Array = probelist;
    probelist = Array.id;
end
np = length(probelist);
S = zeros(1,np);
if ~islogical(selected) && min(selected) > 0
    S(selected) = 1;
elseif max(selected) == 1 && min(selected) == 0
    S(1:length(selected)) = selected;
end
    
j = 1;
while j <= length(varargin)
    if isfield(varargin{j},'X')
        Array = varargin{j};
    elseif strncmp(varargin{j},'fontsize',5)
        j = j+1;
        fontsize = varargin{j};
    elseif strncmp(varargin{j},'single',5)
        singleselection = 1;
    end
    j = j+1;
end

if ~isempty(Array)
    nr = length(unique(Array.X));
    nc = length(unique(Array.Y));
else
    nr = ceil(sqrt(np));
    nc = ceil(np./nr);
end

sz = get(0,'ScreenSize');
F = figure('position',[mean(sz([1 3])) mean(sz([2 4]))+nr.*fontsize.*3 nc.* fontsize.*3 nr.*fontsize.*3]);

p = 0;
for j = 1:nr
    for k = 1:nc
        if isfield(Array,'X')
            p = find(Array.X == j & Array.Y == k);
        else
            p = p+1;
        end
        uicontrol(F,'style','checkbox','string',num2str(probelist(p)),...
            'units','norm',...
            'value',S(p),...
            'position',[(k-1)./nc 1-((j)./(nr+1)) 1./nc 1./(nr+1)],...
            'fontsize',fontsize,...
            'callback',{@HitCheck, p});
    end
end
setappdata(F,'selected',S);
setappdata(F,'single',singleselection);

uicontrol(F,'String','Cancel',...
    'callback',{@HitCheck, 'cancel'},...
    'units','norm',...
    'fontsize',fontsize,...
    'position',[0.01 0.01 0.3 1./(nr+1)]);
if singleselection == 0 
uicontrol(F,'String','Done',...
    'callback',{@HitCheck, 'Done'},...
    'units','norm',...
    'fontsize',fontsize,...
    'position',[0.32 0.01 0.3 1./(nr+1)]);
end
uiwait(F);
selected = getappdata(F,'selected');
close(F);

function HitCheck(a,b, fcn)

F = GetFigure(a);
if strcmp(fcn,'cancel')
    rmappdata(F,'selected');
    uiresume(F);
elseif strcmp(fcn,'Done')
    uiresume(F);
else
    single = getappdata(F,'single');
    selected = getappdata(F,'selected');
    selected(fcn) = get(a,'value');
    setappdata(F,'selected',selected);
    if single
        uiresume(F);
    end
end