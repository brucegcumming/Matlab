%
%bumpdent.m  simple simulation of response to bumps/dents base on zero
%order disparities.
%


% the stimulus wil be 2.0 degrees (+-1) wide. Define the RF over a wider
% range to allow jitter.
position = -1.5:0.05:1.5;

%define response to disparities in range 1:-1
disp = -1:0.05:1;
% disptune is a TI (inverted gaussian) disparity tuning function. In order
% to give it a wide trough at about 20 spikes/sec, as in figure 2b, it is thresholded.
% it is centerd at -0.1 degrees disparity to match the asymmetry in figs 2b
% and 5b

disptune = 60 - 60 .* exp(-((disp+0.1).^2)/(2 * 0.45 .^2));
disptune(find(disptune< 20)) = 20;


%build an RF with disparity tuning defined at each x location
%Three different RFs are build according to the value of simtype
for simtype = [1 2 3];
%
%if simtype == 1, the RF is uniform and larger than the stimulus
%
%if simtype == 2 the RF is uniform but not as wide as the stimulus. At
%values >width and < -width, the RF is set to 0
width = 0.4;

%
%if simtype == 3 each disparity response is weighted as a Gaussian function
%of distance from the RF center
sd = 0.25;

%dispgrid is a matrix defining the response to an combination of position,
%disparity. 

clear dispgrid;

%This for loop builds the three different RF types, depending of simtype
for j = 1:length(position)
    if simtype == 1
        dispgrid(:,j) = disptune';
    elseif simtype == 2
        if abs(position(j)) < width 
            dispgrid(:,j) = disptune';
        else
            dispgrid(:,j) = zeros(length(disptune),1);
        end    
    elseif simtype == 3
        dispgrid(:,j) = disptune' .* exp(-(position(j)^2/(sd^2)));        
    end
    
end

%need to scale the overall values so that the response rates are equivalent
%when the intergral over the RF changes.
scalefactor = mean(disptune)/mean(mean(dispgrid(:,10:end-10)));


%Y = arrary of stimulus disparity values.
%X = array of stimulus position values.
Y = repmat(disp',1,length(position));
X = repmat(position,length(disp),1);

%make pseudocolor plot of RF in  position,disparity space
figure(simtype);
subplot(1,2,1);
hold off;
pcolor(X,Y,dispgrid);
caxis(gca,[0 60]);
shading('interp');
xlabel('Stimulus position');
ylabel('Stimulus disparity');
title('RF response as a function of position,disparity');
hold on;

%bumpstim and dentstim are 1-D stimulus disparity profiles, 1-D functions
%of poistion. The stimulus is 2 degrees wide (cetner of RF and of stimulus at 0)
x = -1:0.05:1;
dentstim = gauss([0 0.4 0.5 -0.25],x);
bumpstim = gauss([0 0.4 -0.5 0.25],x);

ndat = 1;

%jitter can be set between -0.5 and 0.5, otherwise NaN s result
%even with jitter, basic difference between Bumps and Dents preserved.
jitter = 0;

% for each meandisp, caclulate response to Bump, Dent, and plot this
% stimulus over the RF psuedocolor plot
for meandisp = -0.75:0.25:0.75
    meandisps(ndat) = meandisp;
%caluclate response at all values of bumpstim by intepolating dispgrid
    bumpresp(ndat) = mean(interp2(X,Y,dispgrid,x+jitter,bumpstim+meandisp)) .* scalefactor;
    dentresp(ndat) = mean(interp2(X,Y,dispgrid,x+jitter,dentstim+meandisp)) .* scalefactor;
    plot(x,bumpstim + meandisp);
    plot(x,dentstim + meandisp,'r');
    ndat = ndat +1;
end

subplot(1,2,2);
hold off;
%plot bump and dent resplonses, like fig 5b. Superimpose disparity tuning
%like fig 2b
plot(meandisps,bumpresp,'r');
hold on;
plot(meandisps,dentresp,'b');
plot(x,disptune,'k');
yl = get(gca,'Ylim');
set(gca,'Ylim',[0 yl(2)]);
legend('bumps','dents','Normal planar');
xlabel('mean disparity');
ylabel('ratio of spikes to seconds');
end