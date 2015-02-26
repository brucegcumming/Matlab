function Flicker(varargin)

plottype = 'a';
if strcmp(plottype,'overlap')
    steps = 0:100;
    for j = 1:length(steps)
        dx = steps(j);
        G(1,:) = Gauss([-dx 10],[-100:100]);
        G(2,:) = Gauss([dx 10],[-100:100]);
        overlap(j) = sum(min(G));
    end
    hold off;
    plot(overlap);
    hold on;
    plot(1-erf(steps/(10 * sqrt(2))));
else
    
    GetFigure('FlickerChannels');
    logf = 0;
    if logf
    maxf = 6400;
    step = log(maxf)/200;
    sds = exp(0:step:log(maxf));
    else
        maxf = 100;
        sds = 1:maxf/200:maxf
    end
    X = 0:2000;
    for j = 1:length(sds)
        E(j,:) = 1 - Gauss(sds(j),X,'amp',1);
        Fresps(j,:) = Gauss([4000./sds(j) 400./sds(j)],X);
    end
    subplot(2,1,1);
    plot(E');
    xlabel('displacement');
    ylabel('1 - overlap (''Flicker'')');
    subplot(2,1,2);
    plot(E(:,[2 16 32 64 128 256]));    
    xlabel('Filter #');
    ylabel('Flicker');
    title('cross sections');
    GetFigure('Frequency Response');
    dx = 10;
    jumps = [1 2 10 100];
%Fresps is  Filters * jump(tf).    
    for k = 1:length(jumps)
        dx = jumps(k);
        flicker = E(:,dx);  %Flicker for each filter
        %now sum filter responses for each frequency
        for j = 1:size(Fresps,2) %for each stimulus TF
            r = 0;
            for f = 1:size(Fresps,1)
                r(f) = E(f,dx).* Fresps(f,j); %sum response at this TF, scaled by flicker for this fitler
                tr(f) = Fresps(f,j);
            end
            R(k,j) = sum(r);
            Tr(k,j) = sum(tr);
        end
    end
    plot(R');
    set(gca,'xlim',[0 100]);
    max(R);
end