clear all;
close all
load EffectOfThreshold1
ms=10;

for jtype=[EVEN1 ]
    figure('pos',[374   282   905   652])
    for jth=[1 4]
        if jth==1
            subplot(2,2,1)
        elseif jth==4
            subplot(2,2,2)
        end
        % DTC
        curve = squeeze(complex(:,jth,CORR,jtype)) ;
        uncorr = complex(1,jth,UNCORR,jtype);
        
        % Take FT of DTC - UNCORR
        FT = zeros(1,nk);
        for jk=1:nk
            FT(jk) =abs(trapz(disparity,(curve-uncorr)' .* exp(i.*k(jk).*disparity)));
        end
        
        mx=max(curve);
        mn=min(curve);
        amp = max(abs(mx-uncorr),abs(uncorr-mn));
        plot(disparity,(curve+amp-uncorr)/amp,'o','color',col{jth},'MarkerFaceColor',col{jth},'markersize',ms);
        hold on
        curve2 = spline(disparity,curve,disparity2);
        plot(disparity2,(curve2+amp-uncorr)/amp,'color',col{jth},'linewidth',4,'markersize',ms)
        curve = squeeze(complex(:,jth,ANTI,jtype)) ;
        plot(disparity,(curve+amp-uncorr)/amp,'o','color',col{jth},'markersize',ms);
        curve2 = spline(disparity,curve,disparity2);
        plot(disparity2,(curve2+amp-uncorr)/amp,':','color',col{jth})
        axis tight
        xlabel('disparity','fontsize',16)
        set(gca,'ylim',[0 2])
        if jth==1
            title('ODF model','fontsize',16)
            ylabel('Disparity tuning curve','fontsize',16)
        else
            title('New model with high threshold','fontsize',16)
        end
        
        
        % FT:
        if jth==1
            subplot(2,2,3)
        elseif jth==4
            subplot(2,2,4)
        end       
        plot(k,FT,'col',col{jth},'linewidth',4)
        xlabel('frequency','fontsize',16)
        hold on
        axis tight
        lm=get(gca,'ylim');set(gca,'ylim',[lm(1) lm(2)*1.1])
        plot([1 1]*2*pi/lambda,get(gca,'ylim'),'k')
        [mx,indx] = max(FT);
        plot([1 1]*k(indx),get(gca,'ylim'),'col',col{jth},'linestyle',':')
        % Mark on ODF monoc tuning
        plot(k,monoc.*lm(2)./max(monoc),'g');
        ylabel('Fourier transform','fontsize',16)
        set(gca,'ytick',[])
    end
end

for j=1:4;subplot(2,2,j);set(gca,'fontsize',14);end