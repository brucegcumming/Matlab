function [x,DATA] = GetSpikeVals(DATA, ispk, values, dVdt, type, recalc, pcs)
%[x,DATA] = GetSpikeVals(DATA, ispk, values, dVdt, type, recalc, pcs)
% calculate spike properties like Envergy, Var, Width etc

%Do NOT change the order of these definiitons. Cluster params are saved
%as integers....
SPKENERGY=1;
SPKVARE=2;
SPKMAXRATE = 3;
SPKPREMINRATE = 4;
SPKPEAKTIME = 5;
SPKPREMIN = 6;
SPKPEAK = 7;
SPKVAR = 8;
SPKSYMMETRY = 9;
SPKCENTROID = 10;

SPKMINRATE = 11;
SPKMAXRATEA = 12;
SPKMINRATEA = 13;
SPKMEANRATEA = 14;
SPKMAXRATEB = 15;
SPKMINRATEB = 16;
SPKMEANRATEB = 17;
SPKMEANA = 18;
SPKMEANB = 19;
SPKMAXA = 20;
SPKMAXB = 21;
SPKMIN =22;
SPKMINA = 23;
SPKMINB = 24;
SPKHEIGHT = 25;
SPKENERGYA=26;
SPKENERGYB=27;
SPKPREMAXRATE = 28;
TEMPLATEA = 29;
TEMPLATEB = 30;
TEMPLATEC = 31;
SPKSYMNEG = 32;
SPKPCA1= 33;
SPKPCA2= 34;
SPKVARTESTA= 35;
SPKVARTESTB= 36;
SPKPCA3= 37;
SPKPCA4= 38;
SPKPCA12= 39;
SPKMEAN =40;
SPKMAXACCEL = 41;
SPKMINACCEL = 42;
SPKMEANACCEL = 43;
SPKVAREA = 44;
SPKVAREB = 45;
SPKRMSV = 46; %RMS velocity
SPKADC1=47;
SPKADC2=48;
SPKADC3=49;
SPKADC4=50;
ENERGY1 = 51;
ENERGY2 = 52;
ENERGY3 = 53;
ENERGY4 = 54;
if isnan(values) %% return names. Names must match order of actual variable (above). Order can be anything
    if DATA.subprobes > 1
        addsub = 1;
    else
        addsub = 0;
    end
    x = {'Energy' 'Var/Energy' 'MaxRate' 'PreMinRate' 'PeakTime' 'PreMin' 'Peak' 'Var' 'Symmetry' 'Centroid' 'Minrate' 'Maxrate(A)' ... 
    'Minrate(A)' 'Meanrate(A)' 'Maxrate(B)' 'Minrate(B)' 'Meanrate(B)' 'Mean(A)'...
    'Mean(B)' 'Max(A)' 'Max(B)' 'Min' 'Min(A)' 'Min(B)' 'Height' 'EnergyA' 'EnergyB' 'PreMaxRate' ...
    'Template 1' 'Template 2' 'Template 3' '-Symmetry' 'PCA1' 'PCA2' 'test1' 'test2' 'PCA3' 'PCA4' 'PCA1-PCA2' 'Mean' 'AccelMax' 'AccelMin' 'AccelMean'...
    'Var/sqrt(energy)' 'sqrt(Var/Energy)' 'sqrt(Energy)' 'ADC1' 'ADC2' 'ADC3' 'ADC4'};
   DATA = [SPKENERGY SPKVARE SPKMAXRATE SPKVAR SPKPEAK SPKMIN SPKPREMINRATE SPKPEAKTIME SPKPREMIN ...
       SPKSYMMETRY SPKCENTROID SPKMINRATE ...
       SPKMAXRATEA SPKMINRATEA SPKMEANRATEA SPKMAXRATEB SPKMINRATEB ...
    SPKMEANRATEB SPKMEAN SPKMEANA SPKMEANB SPKMAXA SPKMAXB SPKMIN SPKMINA SPKMINB ...
    SPKHEIGHT SPKENERGYA SPKENERGYB SPKPREMAXRATE TEMPLATEA TEMPLATEB ...
    TEMPLATEC SPKSYMNEG SPKPCA1 SPKPCA2 SPKVARTESTA SPKVARTESTB SPKPCA3 SPKPCA4 SPKPCA12  SPKMEAN SPKMAXACCEL SPKMINACCEL SPKMEANACCEL ...
    SPKVAREA SPKVAREB SPKRMSV SPKADC1 SPKADC2 SPKADC3];
if addsub
DATA = [DATA ENERGY1 ENERGY2  ENERGY3 ENERGY4 ];
x = {x{:} 'Energy 1' 'Energy 2' 'Energy 3' 'Energy 4'};
end
return;
end


p = DATA.probe;
if p > size(DATA.cluster,2)
    DATA.cluster{DATA.currentcluster,p} = {};
end

if DATA.currentcluster > size(DATA.cluster,1) || (~isfield(DATA.cluster{DATA.currentcluster,p},'Arange') & isfield(DATA.cluster{1,p},'Arange'))
    DATA.cluster{DATA.currentcluster,p}.Arange = DATA.cluster{1,p}.Arange;
    DATA.cluster{DATA.currentcluster,p}.Brange = DATA.cluster{1,p}.Brange;
    DATA.cluster{DATA.currentcluster,p}.Erange = DATA.cluster{1,p}.Erange;
end
arange = DATA.clusterArange;
brange = DATA.clusterBrange;
splen = size(values,2).*size(values,3);
erange = intersect(DATA.clusterErange,1:splen-1);    
    
if ismember(type,[SPKMAXACCEL SPKMINACCEL SPKMEANACCEL])
    acc = diff(values,2,2);
end
%dVdt and values are already just teh values for ispk
%ispk is only sent so that the correct entries in DATA.Spikes are 
%filled in
if type == SPKENERGY
    if recalc
        x  = sum(dVdt(:,erange).^2,2);
        DATA.Spikes.energy(ispk)= x;
    else
        x = DATA.Spikes.energy(ispk);
    end
elseif type == SPKRMSV
        x = sqrt(DATA.Spikes.energy(ispk));
elseif type == ENERGY1
    erange = 1:(splen/4)-1;
    x  = sum(dVdt(:,erange).^2,2);
elseif type == ENERGY2
    erange = (splen/4):(splen/2)-2;
    x  = sum(dVdt(:,erange).^2,2);
elseif type == ENERGY3
    erange = (splen/2)-1:3*splen/4-3;
    x  = sum(dVdt(:,erange).^2,2);
elseif type == ENERGY4
    erange = (splen-splen/4)-3:splen-4;
    x  = sum(dVdt(:,erange).^2,2);
elseif type == SPKENERGYA
    [mins, imins] = min(values'); 
    imins(find(imins<6)) = 6;
    imins(find(imins>26)) = 26;
    for j = 1:size(dVdt,1)
        x(j)  = sum(dVdt(j,imins(j)-5:imins(j)+5).^2,2);
    end
        DATA.Spikes.energy(ispk)= x;
elseif type == SPKVARE
    if recalc
        x = var(values(:,erange)');
        x = x./DATA.Spikes.energy(ispk);
        DATA.Spikes.vw(ispk) = x;
    else
        x = DATA.Spikes.vw(ispk);
    end
elseif type == SPKVAREA
    x = var(values(:,erange)');
    x = x./sqrt(DATA.Spikes.energy(ispk));
elseif type == SPKVAREB
    x = var(values(:,erange)');
    x = sqrt(x)./sqrt(DATA.Spikes.energy(ispk));
elseif type == SPKMAXRATE
    x = max(dVdt(:,:)');
elseif type == SPKPREMINRATE %min rate value preceding max rate value;
    [x,t] = max(dVdt(:,:)');
    for j = 1:length(t)
        x(j) = min(dVdt(j,1:t(j)));
    end
elseif type == SPKPREMAXRATE %MAX rate preceding peak V
    [x,t] = max(values(:,:)');
    t(find(t > size(dVdt,2))) = size(dVdt,2);
    for j = 1:length(t)
        x(j) = max(dVdt(j,1:t(j)));
    end
elseif type == SPKMAXRATEA
    x = max(dVdt(:,arange)');
elseif type == SPKMAXRATEB
    x = max(dVdt(:,brange)');
elseif type == SPKMINRATE
    x = min(dVdt(:,:)');
elseif type == SPKMINRATEA
    x = min(dVdt(:,arange)');
elseif type == SPKMINRATEB
    x = min(dVdt(:,brange)');
elseif type == SPKMEANA
    x = mean(values(:,arange)');
elseif type == SPKMEAN
    x = mean(values');
elseif type == SPKMEANB
    x = mean(values(:,brange)');
elseif type == SPKMEANRATEA
    x = mean(dVdt(:,arange)');
elseif type == SPKMEANRATEB
    x = mean(dVdt(:,brange)');
elseif type == SPKMAXA
    x = max(values(:,arange)');
elseif type == SPKPCA1
    x = pcs(:,1);
elseif type == SPKPCA2
    x = pcs(:,2);
elseif type == SPKPCA3
    x = pcs(:,3);
elseif type == SPKPCA4
    x = pcs(:,4);
elseif type == SPKPCA12
    x = pcs(:,1)./pcs(:,2);
elseif type == SPKMAXB
    x = max(values(:,brange)');
elseif type == SPKMIN
    x = min(values(:,:)');
    if DATA.plot.nodc
        x = x - mean(values');
    end
elseif type == SPKMINA
    x = min(values(:,arange)');
elseif type == SPKMINB
    x = min(values(:,brange)');
elseif type == SPKPREMIN
    [x,t] = max(dVdt(:,:)');
    for j = 1:length(t)
        x(j) = min(values(j,1:t(j)));
    end
elseif type == SPKPEAKTIME
    [x,t] = max(dVdt(:,:)');
    zc = diff(sign(dVdt(:,:)')); 
    len = size(zc,1);
    for j = 1:length(t)
        tm = find(zc(:,j) < 0 & [1:len]' >= t(j));
        if isempty(tm)
            x(j) = 0;
        else
            x(j) = tm(1); % first zero crossing after peakrate
        end
    end
elseif type == SPKPEAK
    [x,t] = max(values(:,:)');
elseif type == SPKHEIGHT
    [x,t] = max(values(:,:)');
    [y,t] = min(values(:,:)');
    x = x-y;
elseif type == SPKVAR
    x = var(values(:,:)');
elseif type == SPKSYMMETRY || type == SPKSYMNEG
  len =  size(dVdt,2);
  cn  = round(sum((dVdt(:,:).^2)*[1:len]',2)./sum(dVdt(:,:).^2,2));
  for j = 1:length(ispk)
  z(j) = mean(values(j,max([cn(j)-5 1]):cn(j)),2);
  y(j) = mean(values(j,cn(j):min([cn(j)+5 len])),2);
  end
  if type == SPKSYMNEG
  x = atan2(-y,-z);
  else
  x = atan2(y,z);
  end
elseif type == SPKMAXACCEL
    x = max(acc');
elseif type == SPKMINACCEL
    x = min(acc');
elseif type == SPKMEANACCEL
    x = mean(abs(acc'));
elseif type == SPKCENTROID
  len =  size(dVdt,2);
  x  = sum((dVdt(:,:).^2)*[1:len]',2)./sum(dVdt(:,:).^2,2);
elseif ismember(type, [TEMPLATEA TEMPLATEB TEMPLATEC])
    j = 1+type-TEMPLATEA;
    x = values * DATA.Templates(j,:)';
elseif ismember(type, [SPKADC1 SPKADC2 SPKADC3 SPKADC4])
    j = 1+type-SPKADC1;
    x = values(:,DATA.adcid(j));
end

x = double(reshape(x,1,length(x)));
