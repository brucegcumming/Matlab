function agbfig(varargin)
    fontsize = 12;
    fontname = 'Helvetica';
    %colors=(jet(7));
    %colors=colors(shuffle(1:10),:);
    %colors(1,:,:)=[1 0 0];
    props.defaultTextboxshapeFontSize = fontsize;
    props.defaultTextboxshapeFontName = fontname;
    figure;
    props.defaultfigureColorMap = gray;
    close(gcf);
    props.defaultFigureColor = [1 1 1];
    props.defaultlegendfontname = fontname;
    props.defaultlegendfontsize = fontsize;
    props.defaultPatchMarkerFaceColor = [0 0 0];
    props.defaultAxesFontName = fontname;
    props.defaultAxesFontSize = fontsize;
    props.defaultAxesLineWidth = 1;    
    props.defaultAxesTickDir = 'out';
    props.defaultAxesTickDirMode = 'manual';
    props.defaultAxesTitleFontWeight = 'normal';
    props.defaultAxesTitleFontSizeMultiplier = 1.1;
    props.defaultAxesBox = 'off';
    props.defaultColorbarFontName = fontname;    
    props.defaultColorbarFontSize = 10;   
    props.defaultHistogramFaceColor = [0 0 0];
    props.defaultHgjavacomponentBackgroundColor = [1 1 1];
    props.defaultHistogramNumBins = 50;    
    props.defaultLineLineWidth = 1.5;        
    props.defaultScatterMarker = '.';
    %props.defaultScatterMarkerFaceColor = [0 0 0];
    props.defaultScatterlinewidth = 1;
    props.defaultTextboxshapeFontName = fontname;
    props.defaultTextboxshapeFontSize = fontsize;   
    props.defaultuicontrolbackgroundcolor=[1 1 1];
    props.defaultAxesFontSize=9;
    %props.defaultAxesColorOrder = colors;
    fields=fieldnames(props);
    for f = 1:length(fields)
        try
            set(groot,fields{f},props.(fields{f}))
        catch
            mssg(0,'Could not set groot propert %s!',fields{f},'color',[1 0 0]);
        end
    end
end
                     
 