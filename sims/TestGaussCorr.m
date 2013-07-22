function Test_Gauss_Corr(varargin)

fct=varargin{1};

figure;

switch fct
  case '1D'
    if nargin<2, dx=10;  else dx=varargin{2}; end
    if nargin<3, nx=1e4; else nx=varargin{3}; end
    if nargin<4, sd=[0.05:0.05:3]; else sd=varargin{4}; end
    x=linspace(-dx,dx,nx);
    col='brgcmykbrmckbrgcmykbrmck'; col=[col col];
    for i=1:length(sd)
      y1=normpdf(x,-0.5,sd(i)); y2=normpdf(x,0.5,sd(i));
      scatter(y1,y2,4,col(1+mod(i,20)));
      aux=corrcoef(y1,y2);
      disp(['Sep: ' num2str(1/sd(i)) ' Corr: ' num2str(aux(1,2))]);
      c(i) = aux(1,2);
    end
    axis image;
    figure;
    plot(1./sd,c);
    xlabel('separation');
    ylabel('corr');
end
