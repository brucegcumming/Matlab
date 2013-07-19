function SetUIFont(F, varargin)

j = 1;
while j <= length(varargin)
    if isstruct(varargin{j})
        Font = varargin{j};
    end
    j = j+1;
end

f = {'FontSize' 'FontName' 'FontWeight' 'FontAngle'};
    
for j = 1:length(f)
    if isfield(Font,f{j})
        set(F,['DefaultUIControl' f{j}],Font.(f{j}));
    end
end
c = get(F,'children');
types = get(c,'type');
id = strmatch('uicontrol',types);
for j = 1:length(f)
    if isfield(Font,f{j})
        set(c(id),f{j},Font.(f{j}));
    end
end


function SetUIFontRecursive(F)
for k = 1:length(F)
for j = 1:length(f)
    if isfield(Font,f{j})
        if isfigure(F(k))
            set(F(k),['DefaultUIControl' f{j}],Font.(f{j}));
        else
            set(F(k),f{j},Font.(f{j}));
        end
    end
end
c = get(F(k),'children');
types = get(c,'type');
id = strmatch('uicontrol',types);
SetUIFont(c(id),Font);
end
