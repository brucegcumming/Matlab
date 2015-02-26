function calcfft(file)
Z = dlmread(file, ' ');
FZ = fft2(Z);
FZ(1,1) = 0;
SFZ = fftshift(FZ);
AZ=abs(SFZ);
dlmwrite('temp.fft', AZ, ' ');
dlmwrite('temp.in', Z, ' ');
fprintf(1,'\n',AZ);
fprintf(1,'%f\n',AZ);
fprintf(1,'Done\n')

