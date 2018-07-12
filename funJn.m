% funJn calcualtes and returns the Jn integral (eq. (15) in [1])
% Jn=funJn(Tn,betaConst,n,numbMC)
% Jn is the J_{n,\beta} integral (scalar) value
% T_n is the n-based SINR threshold value (eq. (17) in [1])
% betaConstant is path-loss exponent
% n is integer parameter
% betaConst, n, and numbMC are scalars. T_n can be a vector
% numbMc is number of sample (Sobol) points for  quasi-MC integration
%
% Author: H.P. Keeler, Inria Paris/ENS, 2014
%
% References:
% [1] H.P. Keeler, B. Błaszczyszyn and M. Karray,
% 'SINR-based coverage probability in cellular networks with arbitrary
% shadowing', ISIT, 2013

function Jn=funJn(Tn,betaConst,n,numbMC)

% Calculates Jn with various integration methods
% n =2 and 3 uses quad and dblquad respectively
% n>3 uses quasi Monte Carlo based on Sobol points
% function is called by funProbCov; see Corollary 7 in [1]
if nargin==3
    numbMC=10^3; %set default number of qMC sample points
end

Jn=zeros(size(Tn));
%%% Use  quadrature methods
if n==3
    for k=1:length(Tn)
        fv=@(v1,v2)(1./((v1.*v2)+Tn(k))).*(1./((v1.*(1-v2))+Tn(k))).*(v1.*v2.*(1-v2).*(1-v1)).^(2/betaConst).*v1.^(2/betaConst+1);
        Jn(k)=dblquad(fv,0,1,0,1); %perform double qudarature
    end
elseif n==2
    for k=1:length(Tn)
        
        fv=@(v1)(v1.*(1-v1)).^(2/betaConst)./(v1+Tn(k));
        Jn(k)=quad(fv,0,1); %perform single qudarature
    end
elseif n==1
    Jn=ones(size(Tn)); %return ones since J_1=1;
else
    %%% Use QMC method
    numbMCn=numbMC;
    %Use sobol points (can also use 'halton')
    q = qrandstream('sobol',n-1,'Skip',1,'PointOrder','graycode');
    qRandAllVi=qrand(q,numbMCn);
    
    %create eta_i values
    etaValues=zeros(numbMCn,n);
    etaValues(:,1)=prod((qRandAllVi(:,1:end)),2);
    for i=2:n-1
        etaValues(:,i)=(1-qRandAllVi(:,i-1)).*prod((qRandAllVi(:,i:end)),2);
    end
    etaValues(:,n)=(1-qRandAllVi(:,n-1));
    
    %create/sample nominator  of integral kernel
    % (only need to create it once)
    numProdv_i=ones(numbMCn,1);
    for i=1:n-1
        viRand=qRandAllVi(:,i);
        numProdv_i=(viRand).^(i*(2/betaConst+1)-1).*(1-viRand).^(2/betaConst).*numProdv_i; %numerator
    end
    for k=1:length(Tn)
        %create/sample denominator of integral kernel
        denomProdv_i=ones(numbMCn,1);
        for i=1:n-1
            denomProdv_i=(etaValues(:,i)+Tn(k)).*denomProdv_i; %denominator term
        end
        denomProdv_i=(etaValues(:,n)+Tn(k)).*denomProdv_i;
        %generate integral kernel and estimate integral
        kernelInt=numProdv_i./denomProdv_i; %integral kernerl
        Jn(k)=(1+n*Tn(k))/n*mean(kernelInt); % perform quasi-MC step
    end
    
end
