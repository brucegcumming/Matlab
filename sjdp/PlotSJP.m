function PlotSJP (filename,figure_handle)
%function PlotSJP (filename)
%interprets .mid data format and plots graph
%Userdata field of figure is set to SAVE_NAME of file.


%%%%%%%%%%% LOAD FILE

out_array = '';
tolerance = 0;
filename
file_handle = fopen (filename,'r');
if (file_handle~=-1)
     	this_line = fgetl (file_handle);
	while ((tolerance<7))
		if (this_line~=-1)
			out_array = strvcat (out_array,this_line);
			tolerance = 0;

		else
			tolerance = tolerance+1;
		end

		this_line =fgetl (file_handle);
	end
	fclose (file_handle);
	DATAFILE = out_array;      
else
	errordlg ('Load File Failed','File Problem');
	r = 'Load_Failed';
	return;
end



%%%%%%%%%%%% SET DEFAULTS

if (nargin<2)
	hold_status=0;
else
	hold_status=1;
end


HOLD_Y_FLAG = 0;			%flag to see if x axis is held at 0.
LOGX_FLAG = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% DEFAULT PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_ERRORS	= -99;
TEX_SIZE 	= 10;
PLOT_SIZE	= [11.71 10.25];
MARKER_SET 	= 1;
MARKER_SIZE	= 5;
PLOT_NAME	= -99;
FILE_NAME	= -99;
CELL_ID		= -99;
X_AXIS_LABEL	= -99;
Y_AXIS_LABEL 	= -99;
Z_AXIS_LABEL	= -99;
DIAG_HIST	= -99;
FONT_SIZE 	= 12;
PLOT_TYPE 	= -99;
DATA_NAME	= -99;
HGRAM_START	= -99;
HGRAM_NO_BINS	= -99;
HGRAM_BIN_WIDTH = -99;
LEGEND_TEXT     = -99;
DATA_ALL_MEAN   = -99;
DATA_MEAN	= -99;
PLOT_COMMAND    = -99;
SPLINE_INTERPOLATION = -99;
PLANE 		= -99;
CGAUSS_PLOT     = -99;
CGAUSS2_PLOT     = -99;
SINE_PLOT	= -99;
GABOR5_PLOT	= -99;
GABOR_PLOT	= -99;
GABOR_FIX_POS_PLOT	= -99;
GABOR_FIX_PHA_PLOT	= -99;
GAUSSIAN_PLOT		= -99;
GAUSSIAN_OR_PLOT	= -99;
GAUSSIAN_OR_180_PLOT	= -99;
GAUSSIAN_LOG_PLOT	= -99;
GAUSSIAN_FP_PLOT        = -99;
LOG_GAUSSIAN_STATUS  	= -99;
LMONOC_RESP  	= -99;
RMONOC_RESP	= -99;
BLANK		= -99;
UNCORR		= -99;
RIGHT_AXIS	= -99;
RIGHT_AXIS_LABEL = -99;
ERROR_BARS	= 'myerrors';
SQRT_PLOT       = 'OFF';
RUNNING_XMEAN = -99;
RUNNING_XMEDIAN = -99;
SLIDE_FORMAT = -99;
MidnightBlue=[0.1,0.1,0.4]; 
DeepBlue=[0.01,0.01,0.2]; 
AXIS_LABEL_COLOR = [0.1 0.1 0.1];
TITLE_COLOR = [0.1 0.1 0.1];
CURVE_FIT_STYLES = 'k-';
BackColor = [1 1 1];
Error_Bar_Colours = [0 0 0; 0 0 0; 0 0 0; 0 0 0];
PRINT = -99;
BIN_WIDTH = -99;
LEGEND_FONT_SIZE = 9;
COLOR_ORDER=[1 0 0;1 1 0;0 1 1;1 0 1;0 0 1;0 1 0;0.25 0.25 0.25];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% TIP ALL VARIABLES %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contents = DATAFILE;
[contents_y contents_x] = size (contents);
count = 1;
all_fields = [];
all_index = 1;

contents_position = 1;
		
while (contents_position <= contents_y)
        [firstfield remainder] = strtok (contents (contents_position,:));
       	if ((firstfield (1,1)=='*')|(firstfield(1,1)=='#')|(strcmp(firstfield,'DATA')));
		contents_position  = contents_position+1;		
	else
		Value = (sscanf(remainder,'%f'))';
		dims = size (Value);
		variable_index =1;
		if (dims(1,1)~=0)     %is it a float?
			found_a_float = 1;
			command_string = [char(firstfield) '(' num2str(variable_index) ',1:' num2str(dims(1,2)) ')' '=[' num2str(Value) '];' ] ;
			eval (command_string);
			[field_y field_x] = size (firstfield);
			all_fields (all_index,1:field_x) = (firstfield);
			all_index = all_index+1;
			same_flag = 0;
			contents_position = contents_position+1;
			while ((same_flag == 0) & (contents_position<=contents_y))		
				[nextfield remainder] = strtok (contents(contents_position,:));
				if (strcmp (nextfield,firstfield)~=1)
					same_flag=1;	
					contents_position = contents_position -1;
				else	
					variable_index= variable_index+1;
					Value = (sscanf(remainder, '%f'))';
					dims  = size (Value);
                                        command_string = [char(firstfield) '(' num2str(variable_index) ',1:' num2str(dims(1,2)) ')' '=[' num2str(Value) '];' ];
					eval (command_string);
				end
				contents_position = contents_position+1;
			end				
		else		
			found_a_char =1;
			Value = deblank (remainder);	
			while ((double(Value(1,1))==9)|(double(Value(1,1))==32))
				[valuey valuex] = size (Value);
				Value = Value (1:valuey,2:valuex);
			end	
			if (firstfield(1,1)=='#')
				contents_position = contents_position+1;
			else
				dims = size (Value);
				command_string = ['clear ' char(firstfield) ';'];
				eval (command_string);
				command_string = [char(firstfield) '(' num2str(variable_index) ',1:' num2str(dims(1,2)) ')= Value ' ';'];
				eval (command_string);
				[fieldy fieldx] = size (firstfield);
				all_fields (all_index,1:fieldx) = (firstfield);
				all_index = all_index+1;
				same_flag = 0;
				contents_position = contents_position+1;
				while ((same_flag == 0) & (contents_position<=contents_y))		
					[nextfield remainder] = strtok (contents(contents_position,:));
					if (strcmp (nextfield,firstfield)~=1)
						same_flag=1;	
						contents_position = contents_position -1;
					else	
						variable_index= variable_index+1;
						Value = deblank(remainder);
						while ((double(Value(1,1))==9)|(double(Value(1,1))==32))
							[valuey valuex] = size (Value);
							Value = Value (1:valuey,2:valuex);
						end	
						dims  = size (Value);
						command_string = [char(firstfield) '(' num2str(variable_index) ',1:' num2str(dims(1,2)) ')' '= Value' ';'] ;
						eval (command_string);						
					end
					contents_position = contents_position+1;
				end									
			end				
		end
	end
end

	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% MARKER STYLES %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scatter_marker_styles = strvcat ('s','s','^','^','o','o','d','d','p','p','h','h');
marker_edge_styles = strvcat    ('k','k','k','k','k','k','k','k','k','k','k','k');
marker_face_styles = strvcat    ('k','w','k','w','k','w','k','w','k','w','k','w');
line_styles = strvcat ('bo-','mx-','g*-','r+-','ks-','cd','bh-','mp-','r<-','k>-','c^-');
fit_styles = strvcat ('b','m','g','r','k','y');
time_line_styles = strvcat('b-','r-','g-','m-','k-');

if (MARKER_SET==2)
	scatter_marker_styles = strvcat ('s','s','o','o','^','o');
end



if (MARKER_SET==10)
	scatter_marker_styles = strvcat ('s','^','o','s','^','o');
	marker_edge_styles = strvcat ('k', 'k','k','k','k','k');
	marker_face_styles = strvcat ('w','k','w','k','w','k');
	line_styles = strvcat ('bo-','mo-','go-','ro-','ko-','cd');
	fit_styles = strvcat ('b','m','g','r','k','y');
end
if (MARKER_SET==0)
	scatter_marker_styles = strvcat ('s','^','o','s','^','o');
	marker_edge_styles = strvcat ('k','k','k','k','k','k');
	marker_face_styles = strvcat ('w','k','w','k','w','k');
	line_styles = strvcat ('b-','m-','g-','r-','k-','cd');
	fit_styles = strvcat ('b','m','g','r','k','y');
end
if (MARKER_SET==3)
	scatter_marker_styles = strvcat ('s','^','o','s','^','o');
	marker_edge_styles = strvcat ('k', 'k','k','k','k','k');
	marker_face_styles = strvcat ('k','k','k','w','w','w');
	line_styles = strvcat ('bx-','mo-','g*-','r+-','ks-','cd');
	fit_styles = strvcat ('b','m','g','r','k','y');
end
if (MARKER_SET==42)
	
	scatter_marker_styles = strvcat ('s','s','s','s','s','s','s','s');
	marker_edge_styles = strvcat ('k', 'k','k','k','k','k','k','k');
	marker_face_styles = strvcat ('k','w','k','k','w','w','w','k');
	line_styles = strvcat ('ks-','ks-','ks-','ks-','ks-','ks-','ks-','ks-');
	fit_styles = strvcat ('b','m','g','r','k','y');
end

SAVE_NAME = Midas_Get_Field_No_Rep (DATAFILE,'SAVE_NAME');
if (ischar (SAVE_NAME)~=1)
     return;
end



if (PLOT_NAME==-99) PLOT_NAME = FILE_NAME; end;
if (PLOT_NAME==-99)  PLOT_NAME = CELL_ID;  end;
PLOT_NAME = strrep (PLOT_NAME,'_',' ');
if (strcmp (PLOT_NAME,'none')==1)   PLOT_NAME = ''; end

X_AXIS_LABEL = strrep (X_AXIS_LABEL,'_',' ');
Y_AXIS_LABEL = strrep (Y_AXIS_LABEL,'_',' ');
if (X_AXIS_LABEL==-99)     X_AXIS_LABEL = ''; end;
if (Y_AXIS_LABEL==-99)     Y_AXIS_LABEL = ''; end;

DATA_ARRAY = Midas_Get_Data (DATAFILE, 'DATA');

%DATA_ARRAY(:,5,:) = DATA_ARRAY(:,3,:); %converts from standard errors to standard deviations.


if (strcmp(SQRT_PLOT,'ON')==1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% OPEN FIGURE AND SET DEFAULTS %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





if (hold_status==0)
	figure_handle = figure;
end;
set (figure_handle,'DefaultLineMarkerSize',MARKER_SIZE);
set (figure_handle,'DefaultTextFontUnits','points');
set (figure_handle,'DefaultTextFontSize',FONT_SIZE);
set (figure_handle,'DefaultAxesFontUnits','points');
set (figure_handle,'DefaultAxesFontSize',FONT_SIZE);
%set (figure_handle,'DefaultAxesColorOrder',COLOR_ORDER);


if (SLIDE_FORMAT==1)
       	set(figure_handle,'InvertHardCopy','off');
	set(figure_handle,'Color','k');
	whitebg(MidnightBlue);
	colordef none;
	whitebg(MidnightBlue);	
	FontSize = 24;
	BackColor = MidnightBlue;
	set(figure_handle,'DefaultTextFontName','Times');
	set(figure_handle,'DefaultAxesFontName','Times');
	set(figure_handle,'DefaultTextFontSize',16);
	set(figure_handle,'DefaultTextColor',[1 1 1]);
	set(figure_handle,'DefaultAxesColor',MidnightBlue);
	set(figure_handle,'DefaultLineColor',[0 1 1]);
	set(figure_handle,'DefaultAxesYColor',[1 1 0]);
	set(figure_handle,'DefaultAxesXColor',[1 1 0]);
	set (figure_handle,'DefaultPatchEdgeColor',[0 0 0]);
	PLOT_SIZE = [400 300 400 300];
	Position_Vector = [400 300 400 300];
	AXIS_LABEL_COLOR = [1 1 1];
	TITLE_COLOR = [1 1 1];
	scatter_marker_styles = strvcat ('s','s','^','^','o','o');
	marker_edge_styles = strvcat ('y', 'c','y','c','y','c');
	marker_face_styles = strvcat ('y','c','y','c','y','c');
	line_styles = strvcat ('yo-','wx-','c*-','r+-','ms-','cd');
	fit_styles = strvcat ('y','w','c','r','m','y');
	CURVE_FIT_STYLES = 'c-';
	set (figure_handle,'DefaultLineMarkerSize',MARKER_SIZE);
	Error_Bar_Colours = [1 1 0; 0 1 1; 1 1 0; 0 1 1];
else
	set(figure_handle,'MenuBar','none');
end




if (RIGHT_AXIS==-99)
	set (figure_handle,'DefaultAxesPosition',[0.15 0.175 0.825 0.7]);
end
if (hold_status==1)
	set (gcf,'PaperUnits','normalized')
	set (gcf,'PaperType','A4');
	set (gcf,'PaperPosition',[0 0 1 1]);
end;		


if (DIAG_HIST~=-99)
	master_handle = figure_handle;
	subplot(2,2,3);
	axis square
end

PLOT_SIZE = PLOT_SIZE *1024/30;
Position_Vector = [550 400 PLOT_SIZE(1,1) PLOT_SIZE(1,2)]

if (hold_status~=1)
     set (figure_handle,'Position',Position_Vector);
end

%#########################################################
%################## PLOT SPECIFIC STUFF ##################
%#########################################################


if (ischar(PLOT_TYPE)~=1)
     PLOT_TYPE = 'Lines_2d';	%default plot_type if none specified.
end

if (strcmp (PLOT_TYPE,'DISPARITY_TUNING')==1)
	PLOT_TYPE = 'Lines_2d';
	HOLD_Y_FLAG = 1;
	line_styles = strvcat('o-','s-','<-','>-');	
end

if (strcmp (PLOT_TYPE,'DISPARITY_RAMP')==1)
    	%PLOT_TYPE = 'Lines_2d';	
	HOLD_Y_FLAG =1;
	line_styles = strvcat('b-','g-','m-','k-','r-');
end	
if (strcmp (PLOT_TYPE,'FREQUENCY_TUNING')==1)
	PLOT_TYPE = 'Lines_2d';
	HOLD_Y_FLAG = 1;
	FREQUENCY_FLAG=1;
	LOGX_FLAG=1;
	line_styles = scatter_marker_styles;
end
if (strcmp (PLOT_TYPE,'ORIENTATION_TUNING')==1)
     PLOT_TYPE = 'Lines_2d';
     HOLD_Y_FLAG = 1;
     line_styles = scatter_marker_styles; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IF DISPARITY RAMP THEN FLIP REVERSE DATA AND PUT AT MAX CROSS CORR VALUE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PLOT_TYPE
if (strcmp(PLOT_TYPE,'DISPARITY_RAMP')==1)
	FORWARD= DATA_ARRAY(:,2,1)
	FORWARD = FORWARD';
	BACK = DATA_ARRAY  (:,2,2)'
	BACKR = fliplr (DATA_ARRAY (:,2,2)')
	FILTER_SHAPE = [1 1 2 2 3 3 3 2 2 1 1]
	BACKR = conv (BACKR,FILTER_SHAPE);		%blur 
	FORWARD = conv (FORWARD,FILTER_SHAPE);		%blur

	BACKR = fft2(BACKR);
	FORWARD = fft2(FORWARD);
	CROSS_CORR = BACKR.*conj(FORWARD);
	CROSS_CORR = real (ifft2(CROSS_CORR));
	
	MAXIMUM_CORR = maxposn(CROSS_CORR)
	BACKR = flipud(DATA_ARRAY(:,2,2));
	[dimsy dimsx] = size (BACKR)
	MAXIMUM_CORR = MAXIMUM_CORR*-1;
	if (MAXIMUM_CORR<0) 
		MAXIMUM_CORR = MAXIMUM_CORR+dimsy;
	end
	for (rcount = 1:MAXIMUM_CORR)
		temp = BACKR(1,1);
		BACKR(1:dimsy-1) = BACKR(2:dimsy);	%rotate to max posn
		BACKR(dimsy) = temp;
	end 
	BACKR
	DATA_ARRAY(:,2,2) = BACKR;
	PLOT_TYPE = 'Lines_2d';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (strcmp (PLOT_TYPE, 'PSF')==1)	     
	set (figure_handle,'NumberTitle','off');
	set (figure_handle,'Name',SAVE_NAME);
	errorbar (DATA_ARRAY (:,1,1),DATA_ARRAY (:,2,1),DATA_ARRAY (:,5,1),'-o');
	set (gca,'ylim',[0 1]);		
	hold on;
end
if (strcmp (PLOT_TYPE, 'PSF_SIMPLE')==1)	     
	set (figure_handle,'NumberTitle','off');
	set (figure_handle,'Name',SAVE_NAME);
	
	
	[xdims ydims zdims] = size (DATA_ARRAY)
	for zcount = 1:zdims;
		N_field = DATA_ARRAY (:,4,zcount);
		total_no = sum (N_field>0);
		if (zcount==1)
			plot(DATA_ARRAY (1:total_no,1,1),DATA_ARRAY (1:total_no,2,1),'o');
			hold on;
		end
		if (zcount==2)
			plot(DATA_ARRAY (1:total_no,1,2),DATA_ARRAY (1:total_no,2,2),'rs');
		end					
	end;



	set (gca,'ylim',[0 1]);		
	
end



if (strcmp (PLOT_TYPE, 'BAR_GRAPH')==1)	
	if (hold_status~=1)
		set (figure_handle,'Position',Position_Vector);
		set (figure_handle,'NumberTitle','off');
		set (figure_handle,'Name',SAVE_NAME);
	end;
	bar (DATA_ARRAY(:,1,1),'y');
	hold on;
	if (DATA_NAME~=-99)
		set (gca,'XTickLabel',DATA_NAME);     
     	end
	[total_no x z] = size (DATA_ARRAY);
	if ((sum(DATA_ARRAY(1:total_no,2,1))~=-99)&(sum(DATA_ARRAY(1:total_no,3,1))~=-99))
		error_field1 = DATA_ARRAY(1:total_no,2,1);
		error_field2 = DATA_ARRAY(1:total_no,3,1);
		for error_count=1:total_no
			x = error_count*ones(1,2)
			y = [error_field1(error_count) error_field2(error_count)]
			plot (x,y);
			y = error_field1(error_count)*ones(1,2)
			al = axis
			x(1,1) = x(1,1)-((al(1,2)-al(1,1))/160)
			x(1,2) = x(1,2)+((al(1,2)-al(1,1))/160)
			plot(x,y);
			y = [error_field2(error_count)]*ones(1,2)
			plot(x,y);
		end		
	end	
end



if (strcmp (PLOT_TYPE, 'Histogram_1d_Horz')==1)	
	if (hold_status~=1)
		set (figure_handle,'Position',[550 400 350 250]);
		set (figure_handle,'NumberTitle','off');
		set (figure_handle,'Name',SAVE_NAME);
	end;
	[xdims ydims zdims] = size (DATA_ARRAY);
	HGRAM_VECTOR = 0:HGRAM_NO_BINS-1;
	HGRAM_VECTOR = HGRAM_START+HGRAM_BIN_WIDTH/2 + (HGRAM_VECTOR*HGRAM_BIN_WIDTH)
	for zcount = 1:zdims;
		for ycount = 1:xdims
			if ((DATA_ARRAY(ycount,4,zcount)~=-99)&(DATA_ARRAY(ycount,1,zcount)~=-99))
				HGRAM_DATA (ycount,zcount) = DATA_ARRAY (ycount,1,zcount);
			else
				HGRAM_DATA (ycount,zcount) = NaN;
			end;
		end
	end;
	[n xout] = hist(HGRAM_DATA,HGRAM_VECTOR);   
	barh (xout,n,1,'y');
end;

if (strcmp (PLOT_TYPE, 'Histogram_1d_Horz_Log')==1)	
	if (hold_status~=1)
		set (figure_handle,'Position',[550 400 350 250]);
		set (figure_handle,'NumberTitle','off');
		set (figure_handle,'Name',SAVE_NAME);
	end;
	[xdims ydims zdims] = size (DATA_ARRAY);
	
	HGRAM_VECTOR = 0:HGRAM_NO_BINS-1;
	HGRAM_VECTOR = HGRAM_START+HGRAM_BIN_WIDTH/2 + (HGRAM_VECTOR*HGRAM_BIN_WIDTH)
	DATA_ARRAY(:,1,:) = log10 (DATA_ARRAY(:,1,:));
	for zcount = 1:zdims;
		for ycount = 1:xdims
			if ((DATA_ARRAY(ycount,4,zcount)~=-99)&(DATA_ARRAY(ycount,1,zcount)~=-99))
				HGRAM_DATA (ycount,zcount) = DATA_ARRAY (ycount,1,zcount);
			else
				HGRAM_DATA (ycount,zcount) = NaN;
			end;
		end
	end;

	[n xout] = hist(HGRAM_DATA,HGRAM_VECTOR);   
	barh (xout,n,1,'y');

	old_handle = gca;
	old_position = get (gca,'position');
	set (gca,'YTickLabel',[]);
	set (gca,'YTick',[]);
	old_position (1,1) = old_position(1,1);
	al = axis
	new_handle = axes('position',old_position,'XTickMode','Manual','Xtick',[],'TickDir','out');
	axis ([al(1,1) al(1,2) (10^al(1,3)) (10^al(1,4))]);
	set (new_handle,'YScale','log');
	ticks = get (new_handle,'ytick')
	ticks = num2str (ticks')
	set (new_handle,'YTickLabel',ticks)	
	axes (old_handle);
	label_handle = 	get (gca,'ylabel');
	posn = 	get (label_handle,'position')
	
	temp = axis;
	
	posn (1,1) = posn(1,1)-(0.07*(temp(1,2)-temp(1,1)))
	set (label_handle,'position',posn);
	
	hold on;
		
	if (DATA_ALL_MEAN~=-99)
     		if (PLOT_COMMAND~=-99)
    		 	[ploty plotx] = size (PLOT_COMMAND);
    			 for (plot_count = 1:ploty)
				  PLOT_COMMAND (plot_count,:);
			   	  eval (PLOT_COMMAND (plot_count,:),'pause');
			end
		end;
		[datay datax dataz] = size (DATA_ALL_MEAN);
		axis_data = axis;
		min_value = axis_data(1,2)-(axis_data(1,2)-axis_data(1,1))/40;
		max_value = axis_data(1,2);
		x = min_value:(max_value-min_value)/2:max_value;
		y = log10(DATA_ALL_MEAN (1,1))*(x./x);
		plot (x,y,fit_styles (1,:));
	end	
end





if (strcmp (PLOT_TYPE, 'Histogram_1d')==1)	
	if (hold_status~=1)
		set (figure_handle,'Position',Position_Vector);
		set (figure_handle,'NumberTitle','off');
		set (figure_handle,'Name',SAVE_NAME);
	end;
	[xdims ydims zdims] = size (DATA_ARRAY);
	
	HGRAM_VECTOR = 0:HGRAM_NO_BINS-1;
	HGRAM_VECTOR = HGRAM_START+HGRAM_BIN_WIDTH/2 + (HGRAM_VECTOR*HGRAM_BIN_WIDTH)
	for zcount = 1:zdims;
		N_FIELD = DATA_ARRAY(:,4,zcount);
		total_no = sum (N_FIELD>0);
		for ycount = 1:xdims
			if ((DATA_ARRAY(ycount,4,zcount)~=-99)&(DATA_ARRAY(ycount,1,zcount)~=-99)&(ycount<=total_no))
				HGRAM_DATA (ycount,zcount) = DATA_ARRAY (ycount,1,zcount);
			else
				HGRAM_DATA (ycount,zcount) = NaN;
			end;
		end
	end;
	HGRAM_DATA

	if (zdims==1)
		[n xout] = hist(HGRAM_DATA,HGRAM_VECTOR);   
		bar (xout,n,1,'y');
    	else
		[n xout] = hist (HGRAM_DATA,HGRAM_VECTOR)
		temp = bar (xout,n,'stacked');
		[ytemp xtemp] = size(temp)
		for (tempcount =1:xtemp)
			set (temp(tempcount),'FaceColor',COLOR_ORDER(tempcount,:));
     		end
		
	end;


end;


if (strcmp (PLOT_TYPE, 'TIME_COURSE')==1)	
	if (hold_status~=1)
		set (figure_handle,'Position',Position_Vector);
		set (figure_handle,'NumberTitle','off');
		set (figure_handle,'Name',SAVE_NAME);
	end;
	[xdims ydims zdims] = size (DATA_ARRAY);
	
	NO_OF_BINS = round(2500/BIN_WIDTH);
	HGRAM_VECTOR = 0:NO_OF_BINS-1;
	HGRAM_VECTOR = BIN_WIDTH/2 + (HGRAM_VECTOR*BIN_WIDTH);
	for zcount = 1:zdims;
		N_FIELD = DATA_ARRAY(:,4,zcount);
		total_no = sum (N_FIELD>0);
		for ycount = 1:xdims
			if ((DATA_ARRAY(ycount,4,zcount)~=-99)&(DATA_ARRAY(ycount,1,zcount)~=-99)&(ycount<=total_no))
				HGRAM_DATA (ycount,zcount) = DATA_ARRAY (ycount,1,zcount);
			else
				HGRAM_DATA (ycount,zcount) = NaN;
			end;
		end

		
	end;
		[n xout] = hist (HGRAM_DATA,HGRAM_VECTOR);	
	for zcount = 1:zdims
		plot (xout,n(:,zcount),time_line_styles(zcount,:))
		hold on;
	end;		
		
		


end;



%############# if not a histogram then remove items w/-99 ######################

DATA_ARRAY = DATA_ARRAY .* (DATA_ARRAY~=-99.00);


if (strcmp (PLOT_TYPE, 'Lines_2d')==1)
	set (figure_handle,'NumberTitle','off');
	set (figure_handle,'Name',SAVE_NAME);
	[xdims ydims zdims] = size (DATA_ARRAY)
	for zcount = 1:zdims;
		N_field = DATA_ARRAY (:,4,zcount);
		total_no = sum (N_field>0);
		if (ydims==7)
			N_field = DATA_ARRAY (:,5,zcount);
			total_no = sum (N_field>0);
		end
		
		
		
		info = line_styles(zcount);
		marker_edge_styles(zcount);
		marker_face_styles(zcount);
		temp_handle = plot (DATA_ARRAY (1:total_no,1,zcount),DATA_ARRAY(1:total_no,2,zcount),line_styles(zcount,:), ...
				'MarkerEdgeColor',marker_edge_styles(zcount),'MarkerFaceColor',marker_face_styles(zcount))
		set (temp_handle, 'MarkerEdgeColor',marker_edge_styles(zcount));
		set (temp_handle, 'MarkerFaceColor',marker_face_styles(zcount));		
		hold on;
	end;


	if (strcmp(ERROR_BARS,'myerrors')==1)
		if (LOGX_FLAG==1)
			set (gca,'xscale','log');
		end;	
		for (zcount = 1:zdims)
			if (ydims==5)
				if (sum(DATA_ARRAY(1:total_no,5,zcount))~=0)
					N_field = DATA_ARRAY (:,4,zcount);
					total_no = sum (N_field>0);
					error_field = DATA_ARRAY(1:total_no,5,zcount);
					for error_count=1:total_no
						x = DATA_ARRAY(error_count,1,zcount)*ones(1,2);
						ymid = DATA_ARRAY(error_count,2,zcount);
						y = [ymid-error_field(error_count) ymid+error_field(error_count)]
						plot (x,y,'Color',Error_Bar_Colours(zcount,1:3));
						if (strcmp(get(gca,'xscale'),'log')==1)
							y = ymid-error_field(error_count)*ones(1,2);
							al = axis;
							old_factor = (al(1,2)/al(1,1));
							new_factor = 10^(log10(old_factor)/160);						
							x(1,1) = x(1,1)*new_factor;
							x(1,2) = x(1,2)/new_factor;
							plot(x,y,'Color',Error_Bar_Colours(zcount,1:3));
							y = [ymid+error_field(error_count)]*ones(1,2);
							plot(x,y,'Color',Error_Bar_Colours(zcount,1:3));
						else				
							y = ymid-error_field(error_count)*ones(1,2);
							al = axis;
							x(1,1) = x(1,1)-((al(1,2)-al(1,1))/160);
							x(1,2) = x(1,2)+((al(1,2)-al(1,1))/160);
							plot(x,y,'Color',Error_Bar_Colours(zcount,1:3));
							y = [ymid+error_field(error_count)]*ones(1,2);
							plot(x,y,'Color',Error_Bar_Colours(zcount,1:3));
						end
					end
				end
			end	 	
		end
	else
		for zcount = 1:zdims;
			if (ydims==5)
				
				if (sum(DATA_ARRAY(1:total_no,5,zcount))~=0)
					N_field = DATA_ARRAY (:,4,zcount);
					total_no = sum (N_field>0);		
					temp_handle = errorbar (DATA_ARRAY (1:total_no,1,zcount),DATA_ARRAY (1:total_no,2,zcount), ...
						DATA_ARRAY (1:total_no,5,zcount),line_styles(zcount,:));
					set (temp_handle, 'MarkerEdgeColor',marker_edge_styles(zcount));
					set (temp_handle, 'MarkerFaceColor',marker_face_styles(zcount));	
	
				end;
			end;
			if (ydims==7)
				N_field = DATA_ARRAY (:,5,zcount)	;
				total_no = sum (N_field>0);		
				temp_handle = errorbar (DATA_ARRAY (1:total_no,1,zcount),DATA_ARRAY (1:total_no,2,zcount,:), ...
				DATA_ARRAY (1:total_no,6,zcount),DATA_ARRAY (1:total_no,7,zcount),line_styles(zcount));
				set (temp_handle, 'MarkerEdgeColor',marker_edge_styles(zcount));
				set (temp_handle, 'MarkerFaceColor',marker_face_styles(zcount));				
			end;
			hold on;
		end;	
	end;

	if (HOLD_Y_FLAG==1)	
		axis_data = axis;		
		axis_data (1,3) = 0;	
		axis (axis_data);	
	end;
	
	       
	if (LOGX_FLAG==1)
		set (gca,'xscale','log');
	end;	
end

if (strcmp (PLOT_TYPE, 'POLAR')==1)	
	set (figure_handle,'NumberTitle','off');
	set (figure_handle,'Name',SAVE_NAME);
	[xdims ydims zdims] = size (DATA_ARRAY)
	for zcount = 1:zdims;
		N_field = DATA_ARRAY (:,4,zcount);
		total_no = sum (N_field>0);
		info = line_styles(zcount);
		marker_edge_styles(zcount);
		marker_face_styles(zcount);
		temp_handle = polar (DATA_ARRAY (1:total_no,1,zcount),DATA_ARRAY(1:total_no,2,zcount),'k');
		hold on;
	end;

end

if (strcmp (PLOT_TYPE, 'PDF')==1)	
	set (figure_handle,'NumberTitle','off');
	set (figure_handle,'Name',SAVE_NAME);
	[xdims ydims zdims] = size (DATA_ARRAY)
	for zcount = 1:zdims;
		N_field = DATA_ARRAY (:,4,zcount);
		total_no = sum (N_field>0);
		info = line_styles(zcount);
		marker_edge_styles(zcount);
		marker_face_styles(zcount);
		temp_handle = plot (DATA_ARRAY (1:total_no,1,zcount),DATA_ARRAY(1:total_no,2,zcount),'k');
		hold on;
	end;


	if (HOLD_Y_FLAG==1)	
		axis_data = axis;		
		axis_data (1,3) = 0;	
		axis (axis_data);	
	end;
		       
	if (LOGX_FLAG==1)
		set (gca,'xscale','log');
	end;	
end







if (strcmp (PLOT_TYPE, 'ScatterPlot_2d')==1)
	if (DIAG_HIST~=-99)
      		set (figure_handle,'Position',[500 300 400 400]);
	end
	   
	set (figure_handle,'NumberTitle','off');
	set (figure_handle,'Name',SAVE_NAME);
	[xdims ydims zdims] = size (DATA_ARRAY)
	
	for zcount = 1:zdims;
		N_field = DATA_ARRAY (:,4,zcount)
		total_no = sum (N_field>0);
		if (ydims==7)

			N_field = DATA_ARRAY (:,5,zcount);
			total_no = sum (N_field>0);
		end
		temp_handle = plot (DATA_ARRAY (1:total_no,1,zcount),DATA_ARRAY(1:total_no,2,zcount),scatter_marker_styles(zcount), ...
				'MarkerEdgeColor',marker_edge_styles(zcount),'MarkerFaceColor',marker_face_styles(zcount))
		
		
		hold on;
	end;




	if (PLOT_COMMAND~=-99)
	     	[ploty plotx] = size (PLOT_COMMAND);	
		for (plot_count = 1:ploty)
		  	PLOT_COMMAND (plot_count,:);
   			eval (PLOT_COMMAND (plot_count,:),'pause');
   		end
	end
	
	al = axis;			
	for zcount = 1:zdims;

		if ((ydims==5)&(strcmp(ERROR_BARS,'myerrors')~=1))
			N_field = DATA_ARRAY (:,4,zcount);
			total_no = sum (N_field>0);		
			temp_handle = errorbar (DATA_ARRAY (1:total_no,1,zcount),DATA_ARRAY (1:total_no,2,zcount), ...
					DATA_ARRAY (1:total_no,5,zcount),scatter_marker_styles(zcount));
			set (temp_handle, 'MarkerEdgeColor',marker_edge_styles(zcount));
			set (temp_handle, 'MarkerFaceColor',marker_face_styles(zcount));	
		end;

		if ((ydims==5)&(strcmp(ERROR_BARS,'myerrors')==1))	
			N_field = DATA_ARRAY (:,4,zcount);
			total_no = sum (N_field>0);		
			Lower_Errs = DATA_ARRAY(1:total_no,2,zcount)-DATA_ARRAY(1:total_no,5,zcount);
			Upper_Errs = DATA_ARRAY(1:total_no,2,zcount)+DATA_ARRAY(1:total_no,5,zcount);
			Mean_X = DATA_ARRAY(1:total_no,1,zcount);
			Mean_Y = DATA_ARRAY(1:total_no,2,zcount);
			
			for (error_count = 1:total_no)
				if (DATA_ARRAY(error_count,5,zcount)~=0)
					x = Mean_X(error_count)*ones(1,2);
					y = [Lower_Errs(error_count) Upper_Errs(error_count)]
					plot (x,y,'Color',Error_Bar_Colours(zcount,1:3));
				
					if (strcmp(get(gca,'xscale'),'log')==1)
							y = Lower_Errs(error_count)*ones(1,2)
							old_factor = (al(1,2)/al(1,1))
							new_factor = 10^(log10(old_factor)/160);    			
							x(1,1) = x(1,1)*new_factor
							x(1,2) = x(1,2)/new_factor
							plot(x,y,'Color',Error_Bar_Colours(zcount,1:3));
							y = Upper_Errs(error_count)*ones(1,2);
							plot(x,y,'Color',Error_Bar_Colours(zcount,1:3));
					else
							y = Lower_Errs(error_count)*ones(1,2)
							x(1,1) = x(1,1)-((al(1,2)-al(1,1))/160)
							x(1,2) = x(1,2)+((al(1,2)-al(1,1))/160)
							plot(x,y,'Color',Error_Bar_Colours(zcount,1:3));
							y = Upper_Errs(error_count)*ones(1,2)
							plot(x,y,'Color',Error_Bar_Colours(zcount,1:3));
					end
				end;
			end;
		end
		
		if ((ydims==7)&(strcmp(ERROR_BARS,'myerrors')~=1))			
			N_field = DATA_ARRAY (:,5,zcount);
			total_no = sum (N_field>0);		
			Lower_Errs = DATA_ARRAY(1:total_no,2,zcount)-DATA_ARRAY(1:total_no,6,zcount);
			Upper_Errs = DATA_ARRAY(1:total_no,7,zcount)-DATA_ARRAY(1:total_no,2,zcount);
			temp_handle = errorbar (DATA_ARRAY (1:total_no,1,zcount),DATA_ARRAY (1:total_no,2,zcount), ...
				Lower_Errs,Upper_Errs,scatter_marker_styles(zcount));
			set (temp_handle, 'MarkerEdgeColor',marker_edge_styles(zcount));
			set (temp_handle, 'MarkerFaceColor',marker_face_styles(zcount));	
		end;

		if ((ydims==7)&(strcmp(ERROR_BARS,'myerrors')==1))			
			N_field = DATA_ARRAY (:,5,zcount);
			total_no = sum (N_field>0);		
			Lower_Errs = DATA_ARRAY(1:total_no,6,zcount);
			Upper_Errs = DATA_ARRAY(1:total_no,7,zcount);
			Mean_X = DATA_ARRAY(1:total_no,1,zcount);
			Mean_Y = DATA_ARRAY(1:total_no,2,zcount);


			X_ERROR = Midas_Get_Data (DATAFILE,'X_ERROR');	
			if (X_ERROR~=-99)
				LEFT_ERRORS = X_ERROR(1:total_no,1,zcount);
				RIGHT_ERRORS = X_ERROR(1:total_no,2,zcount);
			end		
			for (error_count = 1:total_no)
				x = Mean_X(error_count)*ones(1,2);
				y = [Lower_Errs(error_count) Upper_Errs(error_count)];
				plot (x,y,'Color',Error_Bar_Colours(zcount,1:3));
				

				if (strcmp(get(gca,'xscale'),'log')==1)
						y = Lower_Errs(error_count)*ones(1,2);
						old_factor = (al(1,2)/al(1,1));
						new_factor = 10^(log10(old_factor)/160);						
						x(1,1) = x(1,1)*new_factor;
						x(1,2) = x(1,2)/new_factor;
						plot(x,y,'Color',Error_Bar_Colours(zcount,1:3));
						y = Upper_Errs(error_count)*ones(1,2);
						plot(x,y,'Color',Error_Bar_Colours(zcount,1:3));
				
				else
						y = Lower_Errs(error_count)*ones(1,2);
						x(1,1) = x(1,1)-((al(1,2)-al(1,1))/160);
						x(1,2) = x(1,2)+((al(1,2)-al(1,1))/160);
						plot(x,y,'Color',Error_Bar_Colours(zcount,1:3));
						y = Upper_Errs(error_count)*ones(1,2);
						plot(x,y,'Color',Error_Bar_Colours(zcount,1:3));
				end

				if (X_ERROR~=-99)
					y = Mean_Y(error_count)*ones(1,2);
					x = [LEFT_ERRORS(error_count) RIGHT_ERRORS(error_count)];
					plot (x,y,'Color',Error_Bar_Colours(zcount,1:3));

					if (strcmp(get(gca,'yscale'),'log')==1)
						x = LEFT_ERRORS(error_count)*ones(1,2)
						old_factor = (al(1,4)/al(1,3))
						new_factor = 10^(log10(old_factor)/160);      				
						y(1,1) = y(1,1)*new_factor
						y(1,2) = y(1,2)/new_factor
						plot(x,y,'Color',Error_Bar_Colours(zcount,1:3));
						x = RIGHT_ERRORS(error_count)*ones(1,2);	       			
						plot(x,y,'Color',Error_Bar_Colours(zcount,1:3));	
					else
						x = LEFT_ERRORS(error_count)*ones(1,2)
						y(1,1) = y(1,1)-((al(1,4)-al(1,3))/160)
						y(1,2) = y(1,2)+((al(1,4)-al(1,3))/160)
						plot(x,y,'Color',Error_Bar_Colours(zcount,1:3));
						x = RIGHT_ERRORS(error_count)*ones(1,2)
						plot(x,y,'Color',Error_Bar_Colours(zcount,1:3));       			
					end
				end


			end;
		end;

		
	end;



	for zcount = 1:zdims;
		N_field = DATA_ARRAY (:,4,zcount)
		total_no = sum (N_field>0);
		if (ydims==7)
			N_field = DATA_ARRAY (:,5,zcount);
			total_no = sum (N_field>0);
		end
		temp_handle = plot (DATA_ARRAY (1:total_no,1,zcount),DATA_ARRAY(1:total_no,2,zcount),scatter_marker_styles(zcount), ...
				'MarkerEdgeColor',marker_edge_styles(zcount),'MarkerFaceColor',marker_face_styles(zcount));			
		hold on;
	end;

		hold on;
	if (DIAG_HIST~=-99)
		axis square
	end


end

if (RUNNING_XMEAN(1,1)>0)

	[xdims ydims zdims] = size (DATA_ARRAY);
	for zcount = 1:zdims;
		N_field = DATA_ARRAY (:,4,zcount)
		total_no = sum (N_field>0);
		if (ydims==7)
			N_field = DATA_ARRAY (:,5,zcount);
			total_no = sum (N_field>0);
		end
		ALL_DATAX = [ALL_DATAX; DATA_ARRAY(1:total_no,1,zcount)];
		ALL_DATAY = [ALL_DATAY; DATA_ARRAY(1:total_no,2,zcount)];
	end

	al = axis;
	UPPER_LIM = al(1,4);
	LOWER_LIM = al(1,3);
	LIM_RANGE = al(1,4)-al(1,3);

	MEAN_DATA = [];
	Y_POSN = [];

	if (strcmp(get(gca,'xscale'),'log')==1)
		ALL_DATAX = log10( ALL_DATAX );
	end

	for (run_count = 1:RUNNING_XMEAN(1,2))
		this_centre = LOWER_LIM + run_count*LIM_RANGE/RUNNING_XMEAN(1,2);
		this_lower_lim = this_centre - RUNNING_XMEAN (1,1)/2;
		this_upper_lim = this_centre + RUNNING_XMEAN (1,1)/2;
		MEAN_DATA(run_count) = sum((ALL_DATAY>this_lower_lim).*(ALL_DATAY<this_upper_lim).*ALL_DATAX)/ ...
						sum((ALL_DATAY>this_lower_lim).*(ALL_DATAY<this_upper_lim));
		Y_POSN(run_count) = this_centre;
				
	end
	if (strcmp(get(gca,'xscale'),'log')==1)
	 	MEAN_DATA = 10.^MEAN_DATA;
	end
	plot (MEAN_DATA,Y_POSN);
end;


if (RUNNING_XMEDIAN(1,1)>0)
	[xdims ydims zdims] = size (DATA_ARRAY);
	for zcount = 1:zdims;
		N_field = DATA_ARRAY (:,4,zcount);
		total_no = sum (N_field>0);
		if (ydims==7)
			N_field = DATA_ARRAY (:,5,zcount);
			total_no = sum (N_field>0);
		end
		ALL_DATAX = [ALL_DATAX; DATA_ARRAY(1:total_no,1,zcount)];
		ALL_DATAY = [ALL_DATAY; DATA_ARRAY(1:total_no,2,zcount)];
	end

	al = axis;
	UPPER_LIM = al(1,4);
	LOWER_LIM = al(1,3);
	LIM_RANGE = al(1,4)-al(1,3);

	MEAN_DATA = [];
	Y_POSN = [];

	if (strcmp(get(gca,'xscale'),'log')==1)
		ALL_DATAX = log10( ALL_DATAX );
	end

	[ALLY ALLX] = size (ALL_DATAY);

	for (run_count = 1:RUNNING_XMEAN(1,2))
		this_centre = LOWER_LIM + run_count*LIM_RANGE/RUNNING_XMEAN(1,2);
		this_lower_lim = this_centre - RUNNING_XMEAN (1,1)/2;
		this_upper_lim = this_centre + RUNNING_XMEAN (1,1)/2;
		THIS_DATA_SET = (ALL_DATAY>this_lower_lim).*(ALL_DATAY<this_upper_lim).*(ALL_DATAX-10000);
		THIS_DATA_SET = sort(THIS_DATA_SET);
		SET_SIZE = sum(THIS_DATA_SET<-1);
		THIS_DATA_SET = THIS_DATA_SET(1:SET_SIZE);
		THIS_DATA_SET = THIS_DATA_SET +10000;
		

		if (SET_SIZE>0)
			MEDIAN_DATA(run_count) = median (THIS_DATA_SET);
		else 
			MEDIAN_DATA(run_count) = 0;
		end;

		Y_POSN(run_count) = this_centre;				
	end
	if (strcmp(get(gca,'xscale'),'log')==1)
	 	MEDIAN_DATA = 10.^MEDIAN_DATA;
	end
	plot (MEDIAN_DATA,Y_POSN,'r');
end;




if (strcmp (PLOT_TYPE, 'ScatterPlot_3d')==1)
      	set (figure_handle,'Position',Position_Vector);
	set (figure_handle,'NumberTitle','off');
	set (figure_handle,'Name',SAVE_NAME);
	[xdims ydims zdims] = size (DATA_ARRAY)
 	for zcount = 1:zdims;
		N_field = DATA_ARRAY (:,5,zcount);
		total_no = sum (N_field>0);
		temp_handle = plot3 (DATA_ARRAY (1:total_no,1,zcount),DATA_ARRAY(1:total_no,2,zcount),DATA_ARRAY(1:total_no,3,zcount), ...
				scatter_marker_styles(zcount,:),'MarkerEdgeColor','k');
		hold on;
		set (temp_handle, 'MarkerEdgeColor',marker_edge_styles(zcount));
		set (temp_handle, 'MarkerFaceColor',marker_face_styles(zcount));	
	end;

	grid on;
	
	
	ZAXIS_LABEL = strrep (ZAXIS_LABEL,'_',' ');
	zlabel (ZAXIS_LABEL);

	if (PLOT_COMMAND~=-99)
	     [ploty plotx] = size (PLOT_COMMAND);
	     for (plot_count = 1:ploty)
		  PLOT_COMMAND (plot_count,:)	;
	   	  eval (PLOT_COMMAND (plot_count,:),'pause');
	     end
	end

	
	if (SPLINE_INTERPOLATION~=-99)
		Xdivs = SPLINE_INTERPOLATION(1,1);	
		Ydivs = SPLINE_INTERPOLATION(1,2);
		al = axis;
		[newx newy] = 	meshgrid(al(1,1):((al(1,2)-al(1,1))/Xdivs):al(1,2),al(1,3):((al(1,4)-al(1,3))/Ydivs):al(1,4))
		DATA2 = Midas_Get_Field (DATAFILE,'DATA');
		DATA2 = sortrows(DATA2,2);
		x = DATA2(:,2);
		y = DATA2(:,3);
		z = DATA2(:,4);
		spline = griddata(x,y,z,newx,newy,'cubic')
		mesh (newx,newy,spline)
	end
	
	if (PLANE(1,1)~=-99)	
		if (PLOT_COMMAND~=-99)
		     [ploty plotx] = size (PLOT_COMMAND);
		     for (plot_count = 1:ploty)
			  PLOT_COMMAND (plot_count,:)	;
		   	  eval (PLOT_COMMAND (plot_count,:),'pause');
		     end
		end
		al = axis;
		Xdivs = 20;
		Ydivs = 20;
		[newx newy] = 	meshgrid(al(1,1):((al(1,2)-al(1,1))/Xdivs):al(1,2),al(1,3):((al(1,4)-al(1,3))/Ydivs):al(1,4));
		plane_value = (newx.*PLANE(1,1)+newy.*PLANE(1,2))+PLANE(1,3);
		for xcount = 1:Xdivs
			x = newx(1,xcount)*ones(1,2);
			y = [newy(1) newy(Ydivs)];
			
	                z1 = plane_value(1,xcount);
			z2 = plane_value(Xdivs,xcount);
			
			z = [z1 z2];

			plot3(x,y,z);
			hold on;
		end
		for (ycount = 1:Ydivs)
			y = newy(ycount,1)*ones(1,2);
			x = [newx(1,1) newx(Xdivs,Xdivs)];
			z1 = plane_value(ycount,1);
			z2 = plane_value (ycount,Ydivs);
    			z = [z1 z2];

			plot3(x,y,z);
			hold on;
		end
	
	
		for zcount = 1:zdims;
			N_field = DATA_ARRAY (:,5,zcount);
			total_no = sum (N_field>0);
			temp_handle = plot3 (DATA_ARRAY (1:total_no,1,zcount),DATA_ARRAY(1:total_no,2,zcount),DATA_ARRAY(1:total_no,3,zcount), ...
					scatter_marker_styles(zcount,:),'MarkerEdgeColor','k');
			hold on;
			set (temp_handle, 'MarkerEdgeColor',marker_edge_styles(zcount));
			set (temp_handle, 'MarkerFaceColor',marker_face_styles(zcount));	
		end;
	
		clear z;
		for zcount=1:zdims;
			N_field = DATA_ARRAY(:,5,zcount)	
			total_no = sum (N_field>0)
			for (n = 1:total_no)
     				if (DATA_ARRAY(n,3,zcount)~=0)
		     			y = DATA_ARRAY (n,2,zcount)*ones (2,1);
					x = DATA_ARRAY (n,1,zcount)*ones (2,1);
					
					z(1,1) = PLANE(1,1)*x(1,1)+PLANE(1,2)*y(1,1)+PLANE(1,3);
					z(2,1) = DATA_ARRAY(n,3,zcount);
					[x(1,1) y(1,1) z(2,1) z(1,1)]
					plot3 (x,y,z,'-r');	
				end;
			end;
		end;
	

	end
	
	if (PLANE(1,1)==-99)
		for zcount=1:zdims;
			N_field = DATA_ARRAY(:,5,zcount)	
			total_no = sum (N_field>0)
			for (n = 1:total_no)
     				if (DATA_ARRAY(n,3,zcount)~=0)
		     			z  = 0:DATA_ARRAY (n,3,zcount)/2:DATA_ARRAY(n,3,zcount)
					y = DATA_ARRAY (n,2,zcount)*ones (3,1)
					x = DATA_ARRAY (n,1,zcount)*ones (3,1)
					plot3 (x,y,z,'-r');	
				end;
			end;
		end;
	
	end
	

end;

if (hold_status~=1)
	if (ischar(LEGEND_TEXT)==1) 
		[legh objh] = legend (LEGEND_TEXT,0);
		set (findobj(objh,'Type','text'),'FontUnits','Points');
		set (findobj(objh,'Type','text'),'FontSize',LEGEND_FONT_SIZE); 
		if (strcmp (PLOT_TYPE, 'Histogram_1d')==1)
			set(legh,'Visible','Off');
		end
		
		
		%set (legh,'Tag','legend_tag');
		
	end;	
end;



set (figure_handle,'Tag',deblank (SAVE_NAME));
set (gca,'FontSize',FONT_SIZE-1);
title (PLOT_NAME,'FontUnits','points','FontSize',FONT_SIZE,'Color',TITLE_COLOR);
xlabel (X_AXIS_LABEL,'FontUnits','points','FontSize',FONT_SIZE,'Color',AXIS_LABEL_COLOR);
ylabel (Y_AXIS_LABEL,'FontUnits','points','FontSize',FONT_SIZE,'Color',AXIS_LABEL_COLOR);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT FITTED CURVES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (strcmp(CGAUSS_PLOT,'ON')==1)
   	
   	axis_data= axis;
   	min_value = axis_data(1,1);
   	max_value = axis_data(1,2);
   	x = min_value:(max_value-min_value)/50:max_value;
	z_scores = x-CGAUSS_MEAN;
	
   	z_scores= z_scores./CGAUSS_SD;
	dimensions = size (z_scores);
   	for temp_count =1:dimensions(1,2)
   		fitted_probs (1,temp_count) = ZScore_to_Prob(z_scores(1,temp_count));
   	end	
	plot (x,fitted_probs,'r');	
end

if (strcmp(CGAUSS2_PLOT,'ON')==1)   	
   	axis_data= axis;
   	min_value = axis_data(1,1);
   	max_value = axis_data(1,2);
   	x = min_value:(max_value-min_value)/50:max_value;
	z_scores = x-CGAUSS2_MEAN;	
   	z_scores= z_scores./CGAUSS2_SD;
	dimensions = size (z_scores);
   	for temp_count =1:dimensions(1,2)
   		fitted_probs (1,temp_count) = ZScore_to_Prob(z_scores(1,temp_count));
   	end	
	plot (x,fitted_probs,'b--');	
end




if ((DATA_MEAN~=-99)&(strcmp(PLOT_TYPE,'Histogram_1d_Horz_Log')~=1))
     	[datay datax dataz] = size (DATA_MEAN);
	axis_data = axis;
	min_value = axis_data(1,2)-(axis_data(1,2)-axis_data(1,1))/40;
	max_value = axis_data(1,2);
	for (zcount= 1:dataz)
		x = min_value:(max_value-min_value)/2:max_value;
		y = DATA_MEAN (1,1,zcount)*(x./x);
		plot (x,y,fit_styles (zcount,:));
		if (LEGEND_TEXT~=-99)
			text (x(1,1)+(max_value-min_value)*2,y(1,1),LEGEND_TEXT(zcount,:));
		end;
	end
end


if ((DATA_ALL_MEAN~=-99)&(strcmp(PLOT_TYPE,'Histogram_1d_Horz_Log')~=1))
     	[datay datax dataz] = size (DATA_ALL_MEAN);
	axis_data = axis;
	min_value = axis_data(1,2)-(axis_data(1,2)-axis_data(1,1))/40;
	max_value = axis_data(1,2);
	x = min_value:(max_value-min_value)/2:max_value;
	y = DATA_ALL_MEAN (1,1)*(x./x);
	plot (x,y,fit_styles (zcount,:));
%	text (x(1,1)+(max_value-min_value)*2,y(1,1));

end



%##########################################################
%########### PLOT FITTED GABOR IF GABOR_PLOT = ON #######
%##########################################################


if (strcmp(SINE_PLOT,'ON')==1)
   	axis_data= axis;
   	min_value = axis_data(1,1);
   	max_value = axis_data(1,2); 	
	x = min_value:(max_value-min_value)/50:max_value;
	Gabor_Scores=  SINE_YMEAN+SINE_AMP*cos ((2*3.1413*SINE_FREQ*x)+SINE_PHASE);
	plot (x,Gabor_Scores,CURVE_FIT_STYLES(1,:));	
end



if (strcmp(GABOR_PLOT,'ON')==1)
   	axis_data= axis;
   	min_value = axis_data(1,1);
   	max_value = axis_data(1,2); 	
	x = min_value:(max_value-min_value)/50:max_value;
	Gabor_Scores=  Gabor_Trunc (x,GABOR_XMEAN,GABOR_SD,GABOR_FREQ,GABOR_PHASE,GABOR_YMEAN,GABOR_AMP);	
	plot (x,Gabor_Scores,CURVE_FIT_STYLES(1,:));	
end

if (strcmp(GABOR5_PLOT,'ON')==1)
   	axis_data= axis;
   	min_value = axis_data(1,1);
   	max_value = axis_data(1,2);
	posn =1;  	
	x = min_value:(max_value-min_value)/50:max_value;
	Gabor_Scores=  Gabor_Trunc (x,GABOR5_XMEAN,GABOR5_SD,GABOR5_FREQ,GABOR5_PHASE,GABOR5_YMEAN,GABOR5_AMP);	
	plot (x,Gabor_Scores,CURVE_FIT_STYLES(1,:));	
end

if (strcmp(GABOR_FIX_POS_PLOT,'ON')==1)
   	axis_data= axis;
   	min_value = axis_data(1,1);
   	max_value = axis_data(1,2);
	posn =1;  	
	x = min_value:(max_value-min_value)/50:max_value;
	Gabor_Scores=  Gabor_Trunc (x,GABOR_FIX_POS_XMEAN,GABOR_FIX_POS_SD,GABOR_FIX_POS_FREQ,GABOR_FIX_POS_PHASE,GABOR_FIX_POS_YMEAN,GABOR_FIX_POS_AMP);	
	plot (x,Gabor_Scores,CURVE_FIT_STYLES(1,:));	
end

if (strcmp(GABOR_FIX_PHA_PLOT,'ON')==1)
       	axis_data= axis;
   	min_value = axis_data(1,1);
   	max_value = axis_data(1,2);
	posn =1;  	
	x = min_value:(max_value-min_value)/50:max_value;
	Gabor_Scores=  Gabor_Trunc (x,GABOR_FIX_PHA_XMEAN,GABOR_FIX_PHA_SD,GABOR_FIX_PHA_FREQ,GABOR_FIX_PHA_PHASE,GABOR_FIX_PHA_YMEAN,GABOR_FIX_PHA_AMP);	
	plot (x,Gabor_Scores,CURVE_FIT_STYLES(1,:));	
end



if (strcmp(GAUSSIAN_PLOT,'ON')==1)
   	axis_data= axis;
   	min_value = axis_data(1,1);
   	max_value = axis_data(1,2);
	posn =1;
   	x = min_value:(max_value-min_value)/50:max_value;
	GAUSSIAN_Scores=  Gaussian_Trunc (x,GAUSSIAN_XMEAN,GAUSSIAN_SD,GAUSSIAN_YMEAN,GAUSSIAN_AMP);	
	plot (x,GAUSSIAN_Scores,CURVE_FIT_STYLES(1,:));	
end

if (strcmp(GAUSSIAN_FP_PLOT,'ON')==1)
   	axis_data= axis;
   	min_value = axis_data(1,1);
   	max_value = axis_data(1,2);
	posn =1;
   	x = min_value:(max_value-min_value)/50:max_value;
	GAUSSIAN_Scores=  Gaussian_Trunc (x,0,GAUSSIAN_FP_SD,GAUSSIAN_FP_YMEAN,GAUSSIAN_FP_AMP);	
	plot (x,GAUSSIAN_Scores,CURVE_FIT_STYLES(1,:));	
end





if (strcmp(GAUSSIAN_OR_PLOT,'ON')==1)
   	axis_data= axis;
	min_value = 0;
	max_value = 360;
   	axis_data(1,1) = min_value;
   	axis_data(1,2) = max_value;
	axis (axis_data);
	posn =1;
   	x = min_value:(max_value-min_value)/100:max_value;	
	x2 = -360*( (x-GAUSSIAN_OR_XMEAN)>180);
	x3 = 360* ( (GAUSSIAN_OR_XMEAN-x)>180);
	x3 = x3+x2+x
	GAUSSIAN_Scores=  Gaussian_Trunc (x3,GAUSSIAN_OR_XMEAN,GAUSSIAN_OR_SD,GAUSSIAN_OR_YMEAN,GAUSSIAN_OR_AMP);		
	plot (x,GAUSSIAN_Scores,CURVE_FIT_STYLES(1,:));	
end

if (strcmp(GAUSSIAN_OR_180_PLOT,'ON')==1)
       	axis_data= axis;
	min_value = 0;
	max_value = 180;
   	axis_data(1,1) = min_value;
   	axis_data(1,2) = max_value;
	axis (axis_data);
	posn =1;
   	x = min_value:(max_value-min_value)/200:max_value;	
	x2 = -180*( (x-GAUSSIAN_OR_180_XMEAN)>90);
	x3 = 180* ( (GAUSSIAN_OR_180_XMEAN-x)>90);
	x3 = x3+x2+x
	GAUSSIAN_Scores=  Gaussian_Trunc (x3,GAUSSIAN_OR_180_XMEAN,GAUSSIAN_OR_180_SD,GAUSSIAN_OR_180_YMEAN,GAUSSIAN_OR_180_AMP);		
	plot (x,GAUSSIAN_Scores,CURVE_FIT_STYLES(1,:));	
end

if (strcmp(GAUSSIAN_OR_180_PLOT,'ZERO')==1)
       	axis_data= axis;
	min_value = -90;
	max_value = 90;
   	axis_data(1,1) = min_value;
   	axis_data(1,2) = max_value;
	axis (axis_data);
	posn =1;
   	x = min_value:(max_value-min_value)/200:max_value;	
	x2 = -180*( (x-GAUSSIAN_OR_180_XMEAN)>90);
	x3 = 180* ( (GAUSSIAN_OR_180_XMEAN-x)>90);
	x3 = x3+x2+x
	GAUSSIAN_Scores=  Gaussian_Trunc (x3,GAUSSIAN_OR_180_XMEAN,GAUSSIAN_OR_180_SD,GAUSSIAN_OR_180_YMEAN,GAUSSIAN_OR_180_AMP);		
	plot (x,GAUSSIAN_Scores,CURVE_FIT_STYLES(1,:));	
end



if (strcmp(GAUSSIAN_LOG_PLOT,'ON')==1)         
   	axis_data= axis;
   	min_value = axis_data(1,1);
   	max_value = axis_data(1,2);
	posn =1;
   	
	x = min_value:(max_value-min_value)/10000:max_value;
	
	LOG_GAUSSIAN_Scores=  Gaussian_Trunc (log10(x),GAUSSIAN_LOG_XMEAN,GAUSSIAN_LOG_SD,GAUSSIAN_LOG_YMEAN,GAUSSIAN_LOG_AMP);	
	plot (x,LOG_GAUSSIAN_Scores,CURVE_FIT_STYLES(1,:));	
	if ((strcmp (LOG_GAUSSIAN_STATUS,'LOW PASS')==1)|(strcmp(LOG_GAUSSIAN_STATUS,'HIGH PASS')==1))
		temp = axis;
		xposn = temp(1,1)+(0.05*(temp(1,2)-temp(1,1)))
		yposn = temp(1,4)-(0.05*(temp(1,4)-temp(1,3)))
		text (xposn,yposn,LOG_GAUSSIAN_STATUS,'FontSize',6);
	end
end



hold on;

temp = axis;
textposn = temp(1,2)-(0.09*(temp(1,2)-temp(1,1)));
datax = [(temp(1,2)-(0.03*(temp(1,2)-temp(1,1)))) temp(1,2)];


if (LMONOC_RESP~=-99)
	plot (datax,(LMONOC_RESP(1,1).*(datax./datax)),'-');
	text (textposn,LMONOC_RESP(1,1),'L','FontSize',8);
end

if (RMONOC_RESP~=-99)
	if (((LMONOC_RESP/RMONOC_RESP)>0.95)&(LMONOC_RESP/RMONOC_RESP<1.05))
		textposn = textposn-0.03*(temp(1,2)-temp(1,1));
	end
	plot (datax,(RMONOC_RESP(1,1).*(datax./datax)),'-');
	text (textposn,RMONOC_RESP(1,1),'R','FontSize',8);
end

if (BLANK~=-99)
	plot (datax,(BLANK(1,1).*(datax./datax)),'-');
	text (textposn,BLANK(1,1),'B','FontSize',8);
end

if (UNCORR~=-99)
	plot (datax,(UNCORR(1,1).*(datax./datax)),'-');
	text (textposn,UNCORR(1,1),'U','FontSize',8);
end


%###########################################################
%########### EXECUTE STANDARD MATLAB COMMANDS ##############
%###########################################################

hold on

PLOT_COMMAND
if (PLOT_COMMAND~=-99)
     [ploty plotx] = size (PLOT_COMMAND);
     for (plot_count = 1:ploty)
	  PLOT_COMMAND (plot_count,:);
   	  eval (PLOT_COMMAND (plot_count,:),'pause');
     end
end



set (gca,'box','off');
hold on;


ISLOG = get (gca,'yscale');
if (strcmp(ISLOG,'log')==1)
	ticks = get (gca,'ytick');
	ticks = num2str (ticks');
	set (gca,'yticklabel',ticks);
end;

ISLOG = get (gca,'xscale');
if (strcmp(ISLOG,'log')==1)
	ticks = get (gca,'xtick');
	ticks = num2str (ticks');
	set (gca,'xticklabel',ticks);
end;

if (RIGHT_AXIS~=-99)
	if (PLOT_COMMAND~=-99)
		[ploty plotx] = size (PLOT_COMMAND);
		for (plot_count = 1:ploty)
		    PLOT_COMMAND (plot_count,:);
   		    eval (PLOT_COMMAND (plot_count,:),'pause');
   	 	end	
	end
	
	old_handle = gca;
	old_position = get (gca,'position')
	old_position (1,1) = old_position(1,1)+0.01
	al = axis;
	new_handle = axes('position',old_position,'YAxisLocation','Right','XTickMode','Manual','Xtick',[]);
	axis ([al(1,1) al(1,2) (al(1,3)*RIGHT_AXIS) (al(1,4)*RIGHT_AXIS)]);
	set (new_handle,'YScale',get(old_handle,'Yscale'));
	ticks = get (new_handle,'ytick')
	ticks = num2str (ticks')
	set (new_handle,'YTickLabel',ticks)
	get (new_handle)
	new_handle	
	if (RIGHT_AXIS_LABEL~=-99)
		text_handle = ylabel (RIGHT_AXIS_LABEL);

	end
	axes (old_handle);

	if (PLOT_COMMAND~=-99)
		[ploty plotx] = size (PLOT_COMMAND);
		for (plot_count = 1:ploty)
		  PLOT_COMMAND (plot_count,:);
   		  eval (PLOT_COMMAND (plot_count,:),'pause');
   	 	end	
	end

end

DIAG_HIST

if (DIAG_HIST~=-99) 
 	set (findobj(gcf,'Tag','legend'),'Position',[0.17 0.5 0.07 0.05]);
	if (LEGEND_TEXT~=-99)
		leghandle = legend (LEGEND_TEXT,2);
	end	
	set (leghandle,'Position',[0.17 0.37 0.07 0.05]);
	fonthandle = 	findobj (leghandle,'Type','text')
	
	set (fonthandle,'FontUnits','points');
	set (fonthandle,'FontSize',10);

	set (gcf,'Position',[100 100 600 600]);
	old_handle = gca
	ok_Imhere = 1	
	pos = get (gca,'position')
	posave = (pos(1,1)+pos(1,2))/2;
	pos(1,2) = posave;
	pos(1,1) = posave;

	posave = (pos(1,3)+pos(1,4))/2;	
	pos(1,3) = posave;
	pos(1,4) = posave;
	set (gca,'Position',pos);


	

	right = pos (1,1)+pos(1,3);
	top = pos (1,2)+pos(1,4);

	ISLOG = get (gca,'yscale');
	if (strcmp(ISLOG,'log')==1)
		
		ticks = get (gca,'ytick');


		ticks = num2str (ticks');
		set (gca,'yticklabel',ticks);
	end;

	ISLOG = get (gca,'xscale');
	if (strcmp(ISLOG,'log')==1)		
		ticks = get (gca,'xtick');
		ticks = num2str (ticks');
		set (gca,'xticklabel',ticks);
	end;




	al = axis;
	new_axis_lim1 = (al(1,4)/al(1,1))/10
	new_axis_lim2 = (al(1,3)/al(1,2))*10
	
	new_xsize = sqrt (pos(1,3)^2+pos(1,4)^2)
	new_xsize = pos(1,3)*2;	

	subplot (2,2,2)
	DATA = Midas_Get_Field (DATAFILE,'DATA');
	new_data = DATA(:,2)./DATA(:,3);

	if (strcmp(ISLOG,'log')==1)
		HGRAM_START = log10 (new_axis_lim1)
		HGRAM_NO_BINS = DIAG_HIST
		HGRAM_BIN_WIDTH = (log10 (new_axis_lim2)-log10(new_axis_lim1))/DIAG_HIST;
		new_data = log10(new_data);	
	else
		HGRAM_START = new_axis_lim1;
		HGRAM_NO_BINS = DIAG_HIST;
		HGRAM_BIN_WIDTH = ((new_axis_lim2) - (new_axis_lim1))/DIAG_HIST;
	
	end

	HGRAM_VECTOR = 0:HGRAM_NO_BINS-1;
	HGRAM_VECTOR = HGRAM_START+HGRAM_BIN_WIDTH/2 + (HGRAM_VECTOR*HGRAM_BIN_WIDTH);


	gmean = mean (new_data);
	
	[n xout] = hist(new_data,HGRAM_VECTOR)  
	[tempy tempx] = size (n)
	if (n(tempy)~=0)
		n(1) = 0;
		n(tempx) = 0;
		n = abs(n);
	end

	
	bar_handle = bar (xout,n,1,'y');
	set (bar_handle,'EdgeColor',[0 0 0]);

	temp = get (gca,'position');
	temp (1,3) = new_xsize
	temp (1,4) = new_xsize
	temp (1,1) = right-new_xsize/4;
	temp (1,2) = top-new_xsize/4;

	set (gca,'position',temp);
	get (gca,'position')

	al = axis
	
	al (1,1) = min (xout)	
	al (1,2) = max (xout)
	al (1,4) = al(1,4)*2.5;
	axis (al)
	hold on;
	middleud = (al(1,4)-al(1,3))*6/12;

	yscale = al (1,4)-al(1,3)
	xscale = al (1,2)-al(1,1)

	set (gca,'CameraUpVector',[-xscale yscale 0])
	axis (al)
	hold on;	
	set (gca,'box','off');
	set (gca,'TickDir','Out');

	
	diag_axes_posn = [0.377 0.377 0.429 0.429];
	set (gca,'Position',diag_axes_posn);
	
	text_posnx = al(1,1)- (al(1,2)-al(1,1))/8;
	text_posny = al (1,3)+ (al(1,4)-al(1,3))/10;
	text (text_posnx,text_posny,'Frequency','Rotation',[45]);



	patch ([al(1,1)-(al(1,2)-al(1,1))/100 al(1,1)-(al(1,2)-al(1,1))/100 al(1,2) al(1,2)],[al(1,4) middleud*21/20  middleud*21/20 al(1,4)],BackColor,'EdgeColor',BackColor);
	line ([gmean gmean],[0 middleud*9/9]);
	text_line = sprintf ('Mean %4.3f',10^gmean);
	text (gmean-1.5+0.1,middleud*8/9,text_line,'Rotation',[-45]);
	set (gca,'xtick',[-1 0 1 2]);	
	ticks = get (gca,'xtick')
	ticks = 10.^ticks
	ticks = num2str(ticks')
	set (gca,'xticklabel',ticks);

	ylabels = get (gca,'yticklabel');
	yposns = get (gca,'ytick');
	
	[xsize ysize] = size (yposns)
	for (count = 1:ysize)
		if (yposns(count)>middleud)
			[labely labelx] = size(ylabels);
			ylabels(count,1:labelx) = char (32*ones(1,labelx));
		end
	end
	set (gca,'yticklabel',ylabels);
	set (gcf,'PaperUnits','Inches');
	set (gcf,'PaperPositionMode','auto');
end


if (PRINT~=-99)
	[ysize xsize] = size(PRINT);
	for (count = 1:ysize)
		text_string = sprintf('%4.3f\n',eval(PRINT(count,:)))
	end
	temp = axis;
	xposn = temp(1)+(temp(2)-temp(1))/20;
	yposn = temp(4)-(temp(4)-temp(3))/20;
	temp = text(xposn,yposn,text_string);
	set (temp,'FontSize',8);
end









function r = Midas_Get_Field_No_Rep (DATAFILE,FIELD_NAME);
%function r = Midas_Get_Field_No_Rep (DATAFILE,FIELD_NAME);
%searches through 2d array of characters for string at start
%of line and then returns the following argument.
%doesn't replace _'s with spaces.


dims = size (DATAFILE);
position  =1;

found_flag = 0;

while (position<=dims(1,1))
	[firstfield remainder]= strtok (DATAFILE (position,:));
	if (strcmp (firstfield,FIELD_NAME)==1)
	
		%%%%TRY AND READ INTO VECTOR OF FLOATS

		Value = sscanf(remainder,'%f');
		dims = size (Value);
			
		%%%%IF THAT FAILS READ INTO STRING		

		if (dims(1,1)==0)
			Value = sscanf (remainder,'%s');
		end
		%Value = strrep (Value,'_',' ');
		found_flag = 1;
		break;
	end
	position=position+1; 


end

if (found_flag==0)
	r = -99;
	return;
end

r = Value;


function r = Midas_Get_Data (DATAFILE,FIELD_NAME)
%r = Midas_Get_Data (DATAFILE,FIELD_NAME);
%reads data into a  3d array.
%the third dimension is given by the first character of each line.
%this is usually used to represent different data sets!


dims = size (DATAFILE);
position = 1;
found_flag = 0;
y_posn= ones(1,10);
DATA_ARRAY = zeros(1,1,1);

while (position<=dims(1,1))
	[firstfield remainder] = strtok (DATAFILE (position,:));
	if (strcmp(firstfield,FIELD_NAME)==1)
		found_flag =1;
		Value = sscanf (remainder,'%f');
		z_posn = Value (1,1);
		
		data_size1 = size (Value);
		data_size2 = data_size1 (1,1) -1;
	

		DATA_ARRAY (y_posn( 1,z_posn),1:data_size2,z_posn) = Value (2:data_size2+1,1);
		y_posn (1,z_posn)= y_posn(1,z_posn)+1;

	else 
		if (found_flag==1)
			break;
		end
	end
	position = position +1;
end

if (found_flag==0)
	r = -99;
	return;
end

r = DATA_ARRAY;



function r = Midas_Get_Field (DATAFILE,FIELD_NAME);
%function r = Midas_Get_Field (DATAFILE,FIELD_NAME);
%searches through 2d array of characters for string at start
%of line and then returns the following argument.

dims = size (DATAFILE);
position  =1;

found_flag = 0;

while ((position<=dims(1,1))& (found_flag==0))
	[firstfield remainder]= strtok (deblank(DATAFILE (position,:)));
	if (strcmp (firstfield,FIELD_NAME)==1)
		
		end_flag=0;
		
		%%%% IS IT A FLOAT?
		next_char = double (remainder(1,2));


		%%%% IS A FLOAT.
		if ((next_char>44)&(next_char<58)&(next_char~=47))
			out_array = sscanf (remainder,'%f');
			out_array = out_array';
			position = position+1;
			while ((end_flag==0)&(position<=dims(1,1)))
				[firstfield remainder] = strtok (deblank(DATAFILE(position,:)));
				if (strcmp(firstfield,FIELD_NAME)==1)
					this_line = sscanf (remainder,'%f');
				
					[xsize ysize] = size (out_array);
					out_array = [out_array; this_line'];

				else
					end_flag =1; 
				end
				position = position+1;
			end
			if (isempty(out_array))
				out_array = -99;
			end;

		%%%% IS A STRING

		else
		
			%out_array = sscanf (remainder,'%s');
			out_array = deblank (remainder);
			out_array = char (out_array + ((out_array==9)*23));
			out_array = fliplr(deblank(fliplr(out_array)));
			position = position+1;
			while ((end_flag==0)&(position<=dims(1,1)))
				[firstfield remainder] = strtok (deblank(DATAFILE(position,:)));
				if (strcmp(firstfield,FIELD_NAME)==1)
					%this_line = sscanf (remainder,'%s');

					this_line = deblank(remainder);
					this_line = char (this_line + ((this_line==9)*23));
					this_line = fliplr(deblank(fliplr(this_line)));
					out_array = strvcat (out_array, this_line);
				else
					end_flag =1; 
				end
				position = position+1;
			end
			if (strcmp(deblank(out_array),'Inf')==1)
				out_array = -99;
			end;	
		end

		found_flag = 1;
		break;
	end
	position=position+1; 
end




if (found_flag==0)
	r = -99;
	return;
end


r = out_array;






function r = Gabor_Trunc (xposn,MEAN_POSN,SD,FREQUENCY,PHASE,MEAN_LEVEL,AMPLITUDE)%
%function r = Gabor_Trunc (xposn,MEAN_POSN,SD,FREQUENCY,PHASE,MEAN_LEVEL,AMPLITUDE)
%Plots Gabor function with values truncated at 0.




value = 1;
value =  value .* exp(-1*(MEAN_POSN-xposn).*(MEAN_POSN-xposn)./(2*SD*SD));
value = value .* cos (2*3.1413.*(xposn-MEAN_POSN).*FREQUENCY + PHASE);
value = value.*AMPLITUDE;
value = value+MEAN_LEVEL;

value = value.*(value>0);

r= value;

