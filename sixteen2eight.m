function [eight] = sixteen2eight(sixteen)
for i =1:length(sixteen)
    if (sixteen(i)<256)
        eight(i*2-1) = 0;
        eight(i*2) = uint8(sixteen(i));
    else
        eight(i*2-1) = uint8(floor(double(sixteen(i)) / 256.0));
        eight(i*2)   = uint8(sixteen(i) - (double(eight(i*2-1)) * 256.0));
    end
    if(round(sixteen(i))~=1.0 * double(eight(i*2)) + 256 * double(eight(i*2 -1)))
        disp(['something is wrong ' num2str(i) ' - ' num2str(sixteen(i))]);
    end
end