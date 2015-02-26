function SerialTest()

num = 1;
sp = serial('/dev/tty.USA49Wfd33P2.2','BaudRate',9600);
fopen(sp);
fprintf(sp,'ANSW1\n');
fprintf(sp,'SOR0\n');
fprintf(sp,'NET1\n');
fprintf(sp,'BAUD%d\n',9600);
fprintf(sp,'SP%d\n',50);
fprintf(sp,'LL%d\n',150000);
fprintf(sp,'APL1\n');
fprintf(sp,'%dEN\n',num);
pause(0.01);
fprintf(sp,'%dPOS\n',num);
pause(0.01);
s = fscanf(sp,'%s')
d = sscanf(s,'%d');
fprintf(sp,'%dLR1000\n',num);
pause(0.01);
fprintf(sp,'%dM\n',num);
pause(0.01);
for j = 1:30
fprintf(sp,'%dPOS\n',num);
s = fscanf(sp,'%s');
ts(j) = now;
newd(j) = sscanf(s,'%d');
end
fprintf('End pos %s\n',s);

plot(ts-ts(1),newd);
fprintf(sp,'%dDI\n',num);
pause(0.01);
fclose(sp);
delete(sp);
