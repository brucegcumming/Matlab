function AnovaChi(mode,varargin)

nsim=1000;
if mode == 1
    ndata=3; ngroups=2; npara=1;
    data=randn(nsim,ndata,ngroups);
    x=1:ngroups;
    for i=1:nsim
      y=squeeze(data(i,:,:)); % get simulated data for this particular run
      ym=mean(y); 
      ye=std(y)/sqrt(ndata); 
      %ye=1/sqrt(ndata); % i tried this and it destroys the 1-1 relationship between ANOVA and Chi2
      ye=sqrt(mean(sqr(ye))); % this leads to the graphs below
      line=mean(ym./sqr(ye))/mean(1./sqr(ye)); % weighted average that minimizes Chi2
      chi2=sum(((ym-line)./ye).^2);
      chival(i)=chi2;
      p_chi2(i) =1-cdf('chi2',chi2,ngroups-npara);
      [p_anova(i), details] =anova1(y,[],'off');
      fval(i) = details{2,5};
    end
    figure;
    scatter(p_chi2,p_anova,1);
elseif mode == 2
    
        ndata=10; ngroups=8; npara=2;
    data=randn(nsim,ndata,ngroups);
    x=1:ngroups;
    for i=nsim:-1:1
      y=squeeze(data(i,:,:));
      ym=mean(y); ye=std(y)/sqrt(ndata);
      ye=sqrt(mean(sqr(ye)));
      para=polyfit(x,ym,npara-1);
      line=polyval(para,x);
      chi2=sum(sqr((ym-line)./ye));
      p_chi2(i) =1-cdf('chi2',chi2,ngroups-npara);
      p_anova(i)=anova1(y-ones(ndata,1)*line,[],'off');
    end
    figure;
    scatter(p_chi2,p_anova,1);
end


    function y = sqr(x)
        y = x.^2;