function sd = cv2sd(cv)
% returns SD of a cicular Gaussian that produces a give Circular Variance
if cv > 1
    sd = 0;
else

sd = abs(sqrt(-log(cv) .* 2 .* 28.6462.^2));
end
