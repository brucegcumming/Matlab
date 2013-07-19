function [X,Y,Z] = gratingplot(type)
%
%gratingplot(type)
%
%explore what determines response to mixed contrast cases. 
%type = 1 varies normalization for L eye, and ocularity.
%type = 2 varies binocular normalization for crossed, low contrast
%condtion.
% in these cases, Z is Amplitude(crossed)/Amplitude(low)
%  = ratio(crossed/high)/ratio(low/high)
%
%plot 1 is flat. Llow-Rhigh always equals Llow-Rlow,  if only change
%normalization gain of L eye. 
%plot 2 seems only to depend on 
ngain = [1 1 1 1];
ocu = [1 1];

ny = 1;

if type == 1
for lgain  = [1:0.2:4]
    nx = 1;
for loc  = [0.2:00.2:2]
    norm = [1 1 1];
    ocu(1) = loc;
 %   [h,d,r,A] = grating([ocu(1) ocu(2)], norm);
    
%fixed high contrast in R eye. Low contrast in L eye. 
    norm = [lgain 1  ngain(3)];
    [x,d,r,A] = grating([ocu(1)*0.2 ocu(2)], norm);
    
    norm = [lgain 0.2 * ngain(2)  ngain(3)];
    [l,d,r,A] = grating([ocu(1)*0.2 ocu(2) .* 0.2], norm);
    X(nx,ny) = loc;
    Y(nx,ny) = lgain;
    Z(nx,ny) = x(1)/l(1);
    nx = nx+1;
end
    ny = ny+1;
end
elseif ismember(type,[2 3])
    ngain = [5 5 1];
    for again  = [0.4:0.1:1]
        nx = 1;
        for bgain  = [0.2:0.05:1]
            norm = [1 1 1];
            
            %   [h,d,r,A] = grating([ocu(1) ocu(2)], norm);
            
            %fixed high contrast in R eye. Low contrast in L eye.
            norm = [0.2 * ngain(1) 1  again];
            [x,d,r,A] = grating([ocu(1)*0.2 ocu(2)], norm);
            
            norm = [0.2 * ngain(1) 0.2 * ngain(2)  bgain];
            [l,d,r,A] = grating([ocu(1)*0.2 ocu(2) .* 0.2], norm);
            X(nx,ny) = again;
            Y(nx,ny) = bgain;
            Z(nx,ny) = x(1)/l(1);
            nx = nx+1;
        end
        ny = ny+1;
    end
    if type == 3
        [a,b] = min(abs(Z-1));
        for j = 1:size(Y,2)
            x(j) = X(b(j),j);
            y(j) = Y(b(j),j);
        end
        plot(x,y);
        xlabel('Divisor, H-L');
        ylabel('Divisor, Low');
        return;
    end
end

hold off;

[A,B,C] = fillpmesh(X,Y,Z);
pcolor(A,B,C);
colorbar;
colormap('hot');
