[x2d, y2d] = meshgrid(-63:64,-63:64);
sx = 12;
sy=19;
sf = 0.0313; %cylces per pixel

G = exp(-x2d.^2./2./sx^2) .* exp(-y2d.^2./2./sy^2) ...
        .* cos( 2*pi.*x2d.*sf);
    FT = fft2(G);
    imf = real((ifft2(FT.*FT.*FT)));
    
    
    subplot(2,2,1);
    imagesc(G);
    imc = conv2(G,G);
    subplot(2,2,2),
    imagesc(imc);
    subplot(2,2,3),
    imagesc(imf);
    
   
 x = [-63:64];
 y = exp(-x.^2./2./sx^2)  .* cos( 2*pi.*x.*sf);
 FT = fft(y);
 plot(real(ifft(FT.*FT)))