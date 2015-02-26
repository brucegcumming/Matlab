function ScaleWindow(h, scale, varargin)
% resize a window and its contents

p = get(h,'position');
p(3:4) = p(3:4) .* scale';
set(h,'position',p);
c = get(h,'children');
for j = 1:length(c)
    if strcmp('axes',get(c(j),'type'))
        ScaleWindow(c(j),scale);
        set(c(j),'fontsize',get(c(j),'fontsize').*scale);
    else
        prop = get(c(j));
        if isfield(prop,'Position');
            set(c(j),'position',prop.Position.*scale);
        end
        if isfield(prop,'FontSize');
            set(c(j),'fontsize',prop.FontSize.*scale);
        end
    end
end
    
