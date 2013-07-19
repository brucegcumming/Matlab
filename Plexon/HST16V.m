% HST/16V or HST/32V magnitude, phase, and group delay

% set log frequency range for calculation
logf=[-1:0.01:3];
f=10.^logf;
w=2*pi*f;

% R1 * C sets 0.8 Hz cutoff
R1C=1/(2*pi*0.8);
% R2C sets gain (G-1)
R2C=19*R1C;
% transfer function
H=(1+j*w*R2C./(1+j*w*R1C));

% plot magnitude
figure(1);
semilogx(f,20*log10(abs(H)));
xlabel('Freqeuncy (Hz)');
ylabel('Magnitude (dB)');
hold on;
% draw -3dB line at 0.8 Hz
a=axis;
semilogx([a(1),a(2)],[20*log10(20),20*log10(20)],'k:');
semilogx([a(1),a(2)],[20*log10(20)-3,20*log10(20)-3],'k:');
semilogx([0.8,0.8],[a(3),a(4)],'k:');

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