function ContrastNorm(r, varargin)
%Build contast respn fuciton based on normalization model

coffset = 0.0; %equiavalnt to adding a contrast
c = 0.01:0.01:1;
stimc = c;
cb = 0.0; %additional normalization withour drive
j = 1;
colors = mycolors('white');

while j <= length(varargin)
    if strncmp(varargin{j},'offset',4)
        j = j+1;
        coffset = varargin{j};
    end
    j = j+1;
end
for j = 1:length(coffset)
    c = stimc+coffset(j);
    R =  r(1) .* ( c.^r(2)./(r(3).^r(2) + c.^r(2)+cb.^2));
    R (R < 0) = 0;
    plot(stimc,R,'color',colors{j});
    set(gca,'xscale','log');
    hold on;
end