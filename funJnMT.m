% funJnMT calculates and returns the Jn integral (equations (16) and (17)
% in [2]) used for the factorial moment measure of the STINR process
% defined in [2].
% Jn is a n-variable generalization of the J function originally defined
% in [1] (and calculated in funJn.m); see Remark in [2]
% Jn=funJnMT(Tn,betaConst,n,numbMC,IS)
% Jn is the J_{n,\beta} integral (scalar) value
% T_n is the n-valued STINR threshold value
% betaConstant is path-loss exponent
% n is integer parameter
% betaConst, n, and numbMC are scalars. T_n is a vector of length n
% numbMc is number of sample points for  quasi-MC integration
% IS is for importance sampling step of the q-MC method ie IS=1 turns it on
%
% Author: H.P. Keeler, Inria/ENS, Paris and University of Melbourne, 2018
%
% References
% [1] H.P. Keeler, B. Błaszczyszyn and M. Karray,
% 'SINR-based coverage probability in cellular networks with arbitrary
% shadowing', ISIT, 2013
% [2] B. Błaszczyszyn  and H.P. Keeler, 'Studying the SINR process in
% Poisson networks by using its factorial moment measures', IEEE TOIT,
% 2015


function Jn=funJnMT(Tn,betaConst,n,numbMC,IS)
% Calculates Jn with various integration methods
% n = 2 and 3 uses quad and dblquad respectively
% n>3 uses quasi Monte Carlo based on Sobol points
% function is called by funProbCovMT

if nargin==4
    IS=0; %turn importance sampling off
end

Tn=Tn(:); %converts Tn to a column vector
%%% Use  quadrature methods for n=2 and n=3 cases
if n==3
    fv=@(v1,v2)(1./((v1.*v2)+Tn(1))).*(1./((v1.*(1-v2))+Tn(2))).*(v1.*v2.*(1-v2).*(1-v1)).^(2/betaConst).*v1.^(2/betaConst+1)...
        +(1./((v1.*v2)+Tn(1))).*(1./((v1.*(1-v2))+Tn(3))).*(v1.*v2.*(1-v2).*(1-v1)).^(2/betaConst).*v1.^(2/betaConst+1)...
        +(1./((v1.*v2)+Tn(2))).*(1./((v1.*(1-v2))+Tn(3))).*(v1.*v2.*(1-v2).*(1-v1)).^(2/betaConst).*v1.^(2/betaConst+1);
    Jn=integral2(fv,0,1,0,1); %perform double qudarature
elseif n==2
    fv=@(v1)(v1.*(1-v1)).^(2/betaConst)./(v1+Tn(1))...
        +(v1.*(1-v1)).^(2/betaConst)./(v1+Tn(2));
    Jn=integral(fv,0,1); %perform single qudarature
elseif n==1
    Jn=1; %return a one since J_1=1;
else
    %%% Use QMC method
    numbMCn=numbMC;
    %Use sobol points (can also use 'halton')
    q = qrandstream('sobol',n-1,'Skip',1,'PointOrder','graycode');
    qRandAllVi=qrand(q,numbMCn);
    
    %performance importance sampling step with a simple power-law
    %which approximates the Beta distrubiton for large betaConst
    funIS=1;
    if IS==1
        betaValues=(1:(n-1))*(2/betaConst+1)-1;
        betaMatrix=repmat(betaValues,numbMCn,1);
        qRandAllVi=qRandAllVi.^(1./betaMatrix);
        for i=1:n-1
            funIS=funIS.*betaValues(i).*qRandAllVi(:,i).^(betaValues(i)-1);
        end
    end
    
    %create eta_i values
    etaValues=zeros(numbMCn,n);
    etaValues(:,1)=prod((qRandAllVi(:,1:end)),2);
    for i=2:n-1
        etaValues(:,i)=(1-qRandAllVi(:,i-1)).*prod((qRandAllVi(:,i:end)),2);
    end
    etaValues(:,n)=(1-qRandAllVi(:,n-1));
    
    %create/sample denominator of integral kernel
    denomProdv_i=ones(numbMCn,1);
    for i=1:n-1
        denomProdv_i=(etaValues(:,i)+Tn(i)).*denomProdv_i; %denominator term
    end
    denomProdv_i=(etaValues(:,n)+Tn(n)).*denomProdv_i;
    
    %create/sample nominator  of integral kernel
    % (only need to create it once)
    numProdv_i=ones(numbMCn,1);
    for i=1:n-1
        viRand=qRandAllVi(:,i);
        numProdv_i=(viRand).^(i*(2/betaConst+1)-1).*(1-viRand).^(2/betaConst).*numProdv_i; %numerator
    end
    
    %generate integral kernel and estimate integral
    kernelInt=numProdv_i./denomProdv_i./funIS; %integral kernerl
    Jn=(1+sum(Tn))*mean(kernelInt); % perform quasi-MC step and add it to the previous result
    
end

Jn=Jn/n; %final answer
