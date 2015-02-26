function wta(varargin)

x = -6:6;
G = Gauss(2,x);
plot(x,G);
gamma  = 4;
j = 1;
while j <= length(varargin)
    if strncmpi(varargin{j},'gamma',5)
        j = j+1;
        gamma = varargin{j};
    end
    j = j+1;
end
R(1,:) = (G+G).^gamma;
R(2,:) = (max(G)+G).^gamma - max(G).^gamma;
R(3,:) = (sum(G) -G).^gamma;
R(4,:) = -(R(3,:)-max(R(3,:)));
R = myNormalize(R);
plot(x,R);
legend('Resp','Resp+max','Sum-resp');