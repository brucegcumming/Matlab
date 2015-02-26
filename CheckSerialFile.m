function CheckSerialFile(name, varargin)
%Check for errors in serial file, like failures to set dx in manual expts


check.disparity = 0;
check.reward = 0;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'reward',5)
        check.reward = 1;
    end
    j = j+1;
end


txt = scanlines(name);

if check.disparity
dxe = find(strncmp('dx=',txt,3));
dx = find(strncmp('dx',txt,2));

pid = find(ismember(dx-1,dxe)); %lines where dx follows dx=
for j = 1:length(pid)
    dispa = sscanf(txt{dx(pid(j))-1},'dx=%f');
    dispb = sscanf(txt{dx(pid(j))},'dx%f');
    if dispa ~= dispb
        fprintf('Mismatch on line %d: %.4f vs %.4f\n',dx(pid(j)),dispa,dispb);
    end
end
end

if check.reward
    pid = find(strncmp('#STIMC',txt,6));
    for j = 1:length(pid)
        a = sscanf(txt{pid(j)},'#STIMC %d %d %d %d %d %d %d %f Trw=%f');
        bwrw(j) = a(9);
        binocrw(j) = a(8);
    end
    plot(binocrw,brw);
end