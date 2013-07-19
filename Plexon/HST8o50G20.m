% HST/16V or HST/32V magnitude, phase, and group delay

% set log frequency range for calculation
logf=[-1:0.01:3];
f=10.^logf;
w=2*pi*f;

% Gain = 20x
G=20;
% R = 47Meg, C=0.01uF, 0.33 Hz cutoff
RC=47e6*0.01e-6;
% transfer function
H=G*(j*w*RC./(1+j*w*RC));

% plot magnitude
figure(1);
semilogx(f,20*log10(abs(H)));
xlabel('Freqeuncy (Hz)');
ylabel('Magnitude (dB)');
hold on;
% draw -3dB line at 0.33 Hz
a=axis;
semilogx([a(1),a(2)],[20*log10(20),20*log10(20)],'k:');
semilogx([a(1),a(2)],[20*log10(20)-3,20*log10(20)-3],'k:');
semilogx([0.33,0.33],[a(3),a(4)],'k:');

% plot phase
figure(2);
semilogx(f,(180/pi)*angle(H));
xlabel('Freqeuncy (Hz)');
ylabel('Phase (deg)');
hold on;
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