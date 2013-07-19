% PRA2FP300mag
% Plot the frequency response of a PRA2 board

% Gain in dB
Gain=20*log10(50);

% set frequency range
logf=[-1:0.01:3];
f=10.^logf;
w=2*pi*f;

% 4-pole Butterworth high cut filter response
% FP 300 Hz 
% Linear Technlolgy LTC1563-2 4-pole Butterworth
w0=2*pi*2.56e9/8.45e6;
s=-w0*[cos(pi/8)-i*sin(pi/8);
     cos(pi/8)+i*sin(pi/8);
     cos(3*pi/8)-i*sin(3*pi/8);
     cos(3*pi/8)+i*sin(3*pi/8)];
Hn=(w0^length(s));
Hd=(j*w-s(1));
for ind=2:length(s);
     Hd=Hd.*(j*w-s(ind));
 end
HHC=Hn./Hd;

% highpass filter response
% RC filter 0.497 Hz
R1=68.1e3;
C1=4.7e-6;
H=(j*w*R1*C1)./(1+(j*w*R1*C1));
% two identical RC filters in series (isolated) 
% combination is ~ 0.77 Hz
HLC=H.*H;

% total response
H=HLC.*HHC;

% plot magnitude
figure(1);
semilogx(f,Gain+20*log10(abs(H)),'b')
hold on;
%hold on;
% set y-axis lower limit to 40 dB
%a=axis;
%a(3)=40;
%axis(a);
ylabel('Magnitude');
title('PRA2/FP300-G50 Frequency Response');
grid on;

% plot phase
figure(2);
semilogx(f,(360/(2*pi))*unwrap(angle(H)),'b');
hold on;
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
grid on;

% plot group delay
figure(3);
temp=(H(1:length(H)-1)+H(2:length(H)))./2;
groupdelay=log10(-imag((diff(H)./diff(w))./(temp)));
semilogx(f(1:length(f)-1),log10(-imag((diff(H)./diff(w))./(temp))),'b');
hold on;
xlabel('Frequency (Hz)');
ylabel('Log (base 10) group delay (s)');
grid on;