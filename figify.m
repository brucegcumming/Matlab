function figify(fig, ax, varargin)
% figify(fighandle, axhandle)  makes a matlab figure better for importing to PowerPoint, 
%     by: making lines thicker
%     removing the backgroudn white
%     making the font bigger, and sans serif.
%
%  figify(gcf,gca,'fontname','string'); sets the font style. If the font
%  you name doesn't exist, Matlab does nothing (not even a warning). The
% default font is obvious.
%
%  figify(gcf,gca,'nolabel') removes all the axis tick labels
%
%  figify(gcf,gca,'fontsize',pts); sets the font size
%
%  figify(gcf,gca,'linew',pts); sets the line width

%Aadd listfonts and strcmp to match partial names
fontsize = 20;
lineweight = 2;
fontname = 'Comic Sans MS';
forcefill = 0;
if nargin ==1
    ax = gca;
end
j = 1;
while(j < nargin -1)
    if(strncmpi(varargin{j},'fontsiz',5))
        j = j+1;
        fontsize = varargin{j};
    elseif(strncmpi(varargin{j},'linew',5))
        j = j+1;
        lineweight = varargin{j};
    elseif(strncmpi(varargin{j},'fontname',5))
        j = j+1;
        fontname = varargin{j};
    elseif(strncmpi(varargin{j},'fill',5))
        forcefill = 1;
    elseif(strncmpi(varargin{j},'nolabel',5))
        ylabels = {};
        xlabels = {};
    elseif(strncmpi(varargin{j},'paper',5))
        fontname = 'Arial';
    end
    j = j+1;
end


set(ax,'FontName',fontname,'FontWeight','Bold','LineWidth',lineweight);
for j = 1:length(ax)
    t = get(ax(j),'Title');
    set(t,'FontName',fontname,'FontWeight','Bold');
    x = get(ax(j),'Xlabel');
    set(x,'FontName',fontname,'FontWeight','Bold');
    y = get(ax(j),'Ylabel');
    set(y,'FontName',fontname,'FontWeight','Bold');
    if exist('fontsize')
        set(ax(j),'FontSize',fontsize);
        set(x,'FontSize',fontsize);
        set(y,'FontSize',fontsize);
        set(t,'FontSize',fontsize);
    end
end
if exist('ylabels','var')
    set(ax,'YTickLabel',ylabels);
end
if exist('xlabels','var')
    set(ax,'XTickLabel',xlabels);
end

lines = get(fig,'Children');
for j = 1:length(lines)
    h = get(lines(j));
    if isfield(h,'LineWidth')
        set(lines(j),'LineWidth',lineweight);
    end
end

set(ax,'color','none');

for k = 1:length(ax)
    
    lines = get(ax(k),'Children');
    for j = 1:length(lines)
        f = get(lines(j));
        if isfield(lines(j),'color')
            c = get(lines(j),'color');
        end
        %    set(lines(j),'LineWidth',lineweight,'MarkerFaceColor',c);
        set(lines(j),'LineWidth',lineweight);
        if forcefill
            set(lines(j),'MarkerFaceColor',c);
        end
    end
end
