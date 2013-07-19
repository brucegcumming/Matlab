function [para res] = Fit_Gabor_2_AC_Group(x,y1,y2,varargin)

% function para = Fit_Gabor(x,y,init)
% para = [offset amp1 freq phase1 disp var amp2 phase2]
%           1      2    3    4      5   6    7    8
% init and n_data are optional

DEBUG=0; % 0 for no output, 1 for debugging output (figures)

DISP='off'; % 'final' for some, 'iter' for lots of output of lsqcurvefit (text)


 n_data=1;
 j = 1;
 while j <= nargin-3
     vg = varargin{j};
     if isnumeric(varargin{j})
         ndata = varargin{j};
         j = j+1;
     end
     j = j+1;
 end

     

n=length(x);

x=reshape(x,1,n);
if ischar(y2) % evaluate a fit, at x
    para = Gabor_2_AC(ones(size([x x])));
  para = Gabor_2_AC(y1,x);
  return;
end

y1=reshape(y1,1,n);

y2=reshape(y2,1,n);

sr=[y1 y2];
sr(sr<1)=1; % to prevent too much weight to spikerates<1
sigma=sqrt(sr)./sqrt(n_data); % =mean(sqrt(sr)) / sqrt(n)
sig_min=min(sigma(sigma>0));
sigma(sigma==0)=sig_min;
            
Gabor_2_AC(sigma); % setting sigmas

n=length(x);
nf=pi/(x(2)-x(1)); % Nyquist frequency
mf=n*nf/2;

% initial guess for parameters
if nargin<5 | length(init)~=8
    % baseline guess
    xid = find(abs(x) < 1000); %UC is x = 1000. Don't use that X val
    init(1)=mean([y1(1) y1(2) mean(y1) y1(n-1) y1(n)...
                  y2(1) y2(2) mean(y2) y2(n-1) y2(n)]);
    % amplitude guess
    init(2)=max(abs(y1-init(1)));
    init(7)=max(abs(y2-init(1)));
    % frequency guess
    y11=y1-init(1);
    y22=y2-init(1);
    yf1=fft(y11);
    yf2=fft(y22);
    [dummy idx1]=max(abs(yf1(2:n)));
    [dummy idx2]=max(abs(yf2(2:n)));
    init(3)=mean([idx1 idx2])*2*pi/(max(x(xid))-min(x(xid)));
    init(4)=-angle(yf1(idx1+1));
    init(8)=-angle(yf2(idx2+1));
    % gaussian center guess
    init(5)=mean([mean(x(xid).*abs(y1(xid)-init(1)))/mean(abs(y1(xid)-init(1))) ...
                  mean(x(xid).*abs(y2(xid)-init(1)))/mean(abs(y2(xid)-init(1)))]);
    % gaussian width guess
    init(6)=0.5*mean([mean(x(xid).^2.*abs(y1(xid)-init(1)))/mean(abs(y1(xid)-init(1))) ...
                      mean(x(xid).^2.*abs(y2(xid)-init(1)))/mean(abs(y2(xid)-init(1)))]);
end
    
% initial guess for lower and upper bounds
lb=zeros(1,8);
ub=zeros(1,8);
lb(1)=0.5*init(1);
ub(1)=2.0*init(1);
lb(2)=0.5*init(2);
lb(7)=0.5*init(7);
ub(2)=2.0*init(2);
ub(7)=2.0*init(7);
lb(3)=0;
ub(3)=nf/1.5;
lb(4)=-1000; ub(4)=1000; % in the hope it won't get stuck on an artificial periodic bound
lb(8)=-1000; ub(8)=1000; % in the hope it won't get stuck on an artificial periodic bound
lb(5)=x(1);
ub(5)=x(n);
lb(6)=0.25*(x(2)-x(1))^2;
ub(6)=x(n)-x(1);

if DEBUG 
    Debugging_Info('1',init);
end

% do fit with initially guessed frequency

para1=init;
options=optimset('MaxIter',50,'TolFun',0.02,'Display',DISP);
[para1 res1]=Fit_Gabor_2_AC_sub([4 8],x,y1,y2,sigma,para1,lb,ub,options);

if DEBUG 
    Debugging_Info('2 with original frequency',para1);
end

options=optimset('MaxIter',100,'TolFun',2e-3,'Display',DISP);
[para1 res1]=lsqcurvefit(@Gabor_2_AC,para1,x,[y1 y2]./sigma,lb,ub,options);

if DEBUG 
    Debugging_Info('pre-final with original frequency',para1);
end

% do fit with half the initially guessed frequency

para2=init;
para2(3)=0.5*init(3);
%para2(6)=2*init(6);
options=optimset('MaxIter',50,'TolFun',0.02,'Display',DISP);
[para2 res2]=Fit_Gabor_2_AC_sub([4 8],x,y1,y2,sigma,para2,lb,ub,options);

if DEBUG
    Debugging_Info('2 with half of original frequency',para2);
end

options=optimset('MaxIter',100,'TolFun',2e-3,'Display',DISP);
[para2 res2]=lsqcurvefit(@Gabor_2_AC,para2,x,[y1 y2]./sigma,lb,ub,options);

if res1<res2
    para=para1;
    res=res1;
else
    para=para2;
    res=res2;
end

if DEBUG 
    Debugging_Info('pre-final with half of original frequency',para);
end

options=optimset('MaxIter',1000,'TolFun',1e-6,'Display',DISP);
[para res]=lsqcurvefit(@Gabor_2_AC,para,x,[y1 y2]./sigma,lb,ub,options);

if para(3)>nf/2
    disp(['warning: frequency very high: ' num2str(para(3))]);
    disp(['warning: nyquist frequency  : ' num2str(nf)]);
end

if para(3)<0           % frequency negative
    para(3)=-para(3);
    para(2)=-para(2); % amp1
    para(4)=-para(4); % phase1
    para(7)=-para(7); % amp2
    para(8)=-para(8); % phase2
    disp('corrected negative frequency');
end

if para(2)<0           % first amplitude negative
    para(2)=-para(2);
    para(7)=-para(7);
    para(4)=para(4)+pi;
    para(8)=para(8)+pi;
    disp('corrected negative amplitude');
end

para(4)=mod(para(4),2*pi);
para(8)=mod(para(8),2*pi);

if res1<res1
    if norm(para1-para2)/norm(para2)>0.001
        disp('used lower frequency');
        [res1 para1]
        [res2 para2]
    end
end

if DEBUG 
    Debugging_Info('final',para);
end

% DEBUGGING INFO START -------------------------------------------------------
    function Debugging_Info(tag,para)
        disp(['para_debug : ' num2str(para)]);
        Get_Figure(['debug:' tag]);
        plot(x,y1,'b.-'); plot(x,y2,'b.--');
        Gabor_2_AC(1);
        aux=Gabor_2_AC(para,x);
        plot(x,aux(1:n),'r-'); plot(x,aux(n+1:2*n),'r--');
        Gabor_2_AC(sigma);
        plot([para(5) para(5)],[0 max(aux)],'g-');
    end

% DEBUGGING INFO END -------------------------------------------------------

end

function [para res] = Fit_Gabor_2_AC_sub(fit,x,y1,y2,sigma,init,lb,ub,options)

% function para = Fit_Gabor(x,y,init)
% para = [offset amp1 freq phase1 disp var amp2 phase2]
%           1      2    3    4      5   6    7    8

    Gabor_2_AC(sigma);

    if nargin<6
        lb=[];
    end
    if nargin<7
        lb=[];
    end
    if nargin<8
        options=[];
    end

    fix=setdiff([1:length(init)],fit);

    p_fit=init(fit); lb_fit=lb(fit); ub_fit=ub(fit);
    p_fix=init(fix);

    %[p_fit res]=lsqcurvefit(@(p_fit,x) Gabor_2_AC_sub(p_fit,x,p_fix),p,x,[y1 y2]./sigma,[],[],options);
    [p_fit res]=lsqcurvefit(@Gabor_2_AC_sub,p_fit,x,[y1 y2]./sigma,lb_fit,ub_fit,options);

    para(fit)=p_fit;
    para(fix)=p_fix;

    if para(2)<0
        para(2)=-para(2);
        para(4)=para(4)+pi;
        para(8)=para(8)+pi;
        disp('corrected negative amplitude');
    end

    para(4)=mod(para(4)+pi,2*pi)-pi;
    para(8)=mod(para(8)+pi,2*pi)-pi;

        function y = Gabor_2_AC_sub(p_fit,x)
            para(fit)=p_fit;
            para(fix)=p_fix;
            y=Gabor_2_AC(para,x);
        end

end

function y = Gabor_2_AC(para,x)

% function y = Gabor_2_AC(para,x)

    persistent sigma;

    if nargin==1 % set sigma's
        sigma=para;
        y=-1;
        %disp(['Gabor_2_AC: sigma set to: ' num2str(sigma)]);
    else

        if isempty(sigma)
            sigma=1;
            disp('Gabor_2_AC: sigma set to 1!');
        end

        y1=para(1)+para(2)*cos(para(3)*x-para(4)).*exp(-(x-para(5)).^2/2/para(6));

        y2=para(1)+para(7)*cos(para(3)*x-para(8)).*exp(-(x-para(5)).^2/2/para(6));

        y=[y1 y2];

        y(y<0)=0; % simple rectification

        y=y./sigma;
      end

end    
    
function fig=Get_Figure(tag,option,idx)

    if nargin<2
        option='clear';
    end

    if nargin<3
        idx=1;
    end

    if isempty(findstr(option,'intact'))
        option=[option 'clear'];
    end

    figs=findobj('Tag',tag);

    if isempty(figs) | (~isempty(figs) & length(figs)<idx) | ~isempty(findstr(option,'force'))
        fig=figure;
    else
        fig=figs(idx);
        figure(fig(idx));
    end

    if ~isempty(findstr(option,'clear'))
        clf(fig);
    end

    set(fig,'Tag',tag);
    hold on;

    title(tag);

end