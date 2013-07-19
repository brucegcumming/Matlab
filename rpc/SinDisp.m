function dpamp = SinDisp(amp)
%SinDisp calculated amplitude of disparity selective response
%for a single RCP subunit, as a fucntion of the amplitude of the kernel
% response at that frequency (Fourier amplitude of kernel at F.
% amp is a vector of amplitudes, dpamp returns a vector of disp tuned
% amplitudes

plotresp = 1;
colors = mycolors;
step = pi/100;
x = step:pi/100:2*pi;
dstep = pi/10;
dps = dstep:dstep:2*pi;
th = 0;

for j = 1:length(amp)
    for k = 1:length(dps);
% for each RF, resp is k * sinewave. If rfs have same amplitude
% spectrum, k is the same in both eyes
% dps(k) sums any phase difference between RF and phase difference in
% the stim. We are not interested in phase of resp, just amplitude as a
% function of disparity
        lr = amp(j) .* sin(x);
        rr = amp(j) .* sin(x+dps(k));
% threshold before summation
        rr(find(rr < th)) = th;
        lr(find(lr < th)) = th;
        resps(k) = mean((lr + rr).^2);
    end
    if(plotresp)
        plot(dps,resps,'color',colors{k});
        hold on;
    end
    dpamp(j) = max(resps) - min(resps);
end

