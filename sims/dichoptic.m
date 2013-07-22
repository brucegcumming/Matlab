function dichoptic(varargin)

type = 1;
w = 0.5;
j = 1;
while j <= length(varargin)
    if j == 1 & isnumeric(varargin{j})
        type = varargin{j};
    end
    j = j+1;
end
%try out pure dichoptic motion stimuli.

sf = 0.01;
tf = 0.01;
t = repmat(1:256,256,1);
x = repmat([1:256]',1,256);
if type == 1
R = cos(x .* 2 * pi * sf) .* cos(t .* 2 * pi * tf);
L = sin(x .* 2 * pi * sf) .* sin(t .* 2 * pi * tf);
elseif type == 3
Rc = sin(x .* 2 * pi * sf) .* sin(t .* 2 * pi * tf);
R = cos(x .* 2 * pi * sf + t .* 2 * pi * tf) + w .* Rc;
L = sin(x .* 2 * pi * sf) .* sin(t .* 2 * pi * tf);
else
R = cos(x .* 2 * pi * sf + t .* 2 * pi * tf);
L = sin(x .* 2 * pi * sf) .* sin(t .* 2 * pi * tf);
end
subplot(1,3,1);
imagesc(R);
subplot(1,3,2);
colormap('gray');
imagesc(L);
colormap('gray');
subplot(1,3,3);
imagesc(R+L);
colormap('gray');

