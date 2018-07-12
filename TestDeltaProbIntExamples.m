% Calculates delta probability for 2-base station signal combination
% and second-strongest signal removed under arbitrary shadowing/fading
% in a multi-tier cellular network based on a Poisson
% model outlined in [2] (which is similar to [4]). The calculation 
% methods are based on factorial moment measures and were developed in [2].
% The (quasi Monte Carlo) integration method use shigher order Sobol points
% For SIR (ie noise N=0) results are scale-invariant (ie independent
% of base station desnity lambda)
% Increase the variable numbMC, number of sample points, for more accuracy
% but usually 2^10 is sufficient
% A variation of this code was used to produce Figures 6 to 11 in [2]
% Warning: if epsilonPrime is less than 0.05, potentially too many integral
% terms will be required
%
% Author: H.P. Keeler, Inria/ENS Paris, 2018
%
% References:
% [1] H.P. Keeler, B. Błaszczyszyn and M. Karray,
% 'SINR-based coverage probability in cellular networks with arbitrary 
% shadowing', ISIT, 2013
% [2] B. Błaszczyszyn  and H.P. Keeler, 'Studying the SINR process in
% Poisson networks by using its factorial moment measures', IEEE TOIT,
% 2015
% [3] Błaszczyszyn  and H.P. Keeler, 'Equivalence and comparison of random
% heterogeneous networks', presented at WDN-CN, 2013
% [4] Dhillon, Ganti, Adnrews, Baccelli, 'Modeling and analysis of K-tier
% downlink heterogeneous cellular networks', JSAC, 2012

clear all;close all;clc;
numbMC=2^12; %a base-two number is more efficient for Sobol

%network parameters (most can be ignored when there is no noise ie W=0)
lambda=1; %base station density
betaConst=3; %path-loss exponent
K=1; %path-loss constant
W=1/100; %noise term
%log normal parameters
sigmDb=10;
sigma=sigmDb/10*log(10);
ESTwoBeta=exp(sigma^2*(2-betaConst)/betaConst^2);
a=lambda*pi*ESTwoBeta/K^2; %equation (6) in [1] or (4) in [2]
gammaConst=1; %inteference cancellation factor
x=(W/gammaConst)*a^(-betaConst/2);
%x=0; % x=0 corresponds to no noise


nComb=2; %number of signals combined number
epsilonPrime=0.1; %lower threshold value; less than 0.05 requires too
% many (ie more than 20) integrals

tauNumb=25;
tauPrimeValues=linspace(epsilonPrime,.99,tauNumb);
mu_n=0;
i_max=ceil(1/(gammaConst*epsilonPrime));
i_max_temp=i_max;
probDeltaMCComb=zeros(tauNumb,i_max_temp+1); probDeltaMC_IC=probDeltaMCComb;
ithTermMatrix=probDeltaMCComb;

%integration region is often dided into regions labelled A and B
icIndA=1; icIndB=1; %innitiate indicator functions for integration kernel
for i=0:i_max_temp
    %%% Use QMC integration method
    nMoment=i+nComb;
    
    %if nMoment>1
    numbMCn=numbMC; %scale number of points by dimension
    qNumb=2*nComb+i-1; %number of quasi-random ensembles;n for z; n-1 for v
    q = qrandstream('sobol',qNumb,'Skip',1,'PointOrder','graycode');
    qRandAll=qrand(q,numbMCn);
    
    %create eta_i values (which sum to one)
    etaValues=zeros(numbMCn,nMoment);
    qRandAllVi=qRandAll(:,(nComb+1):end);
    
    %performance importance sampling step with a simple power-law
    %which approximates the Beta distrubiton for large betaConst
    IS=1;
    funIS=1;
    if IS==1
        betaValues=(1:(nMoment-1))*(2/betaConst+1)-1;
        betaMatrix=repmat(betaValues,numbMCn,1);
        qRandAllVi=qRandAllVi.^(1./betaMatrix);
        for j=1:nMoment-1
            funIS=funIS.*betaValues(j).*qRandAllVi(:,j).^(betaValues(j)-1);
        end
    end
    %create beta values
    etaValues(:,1)=prod((qRandAllVi(:,1:end)),2);
    for j=2:nMoment-1
        etaValues(:,j)=(1-qRandAllVi(:,j-1)).*prod((qRandAllVi(:,j:end)),2);
    end
    etaValues(:,nMoment)=(1-qRandAllVi(:,nMoment-1));
    %create/sample numerator  of integral kernel
    RviRand=ones(numbMCn,1);
    %create the product of beta probability dentities
    for j=1:nMoment-1
        viRand=qRandAllVi(:,j);
        RviRand=(viRand).^(j*(2/betaConst+1)-1).*(1-viRand).^(2/betaConst).*RviRand; %numerator
    end
    zPrime=qRandAll(:,(1:nComb));
    for k=1:tauNumb
        tauPrime=tauPrimeValues(k);
        zPrimeVolA=1; zPrimeVolB=1;
        %%%create z values
        zPrimeA=zeros(numbMCn,nMoment);
        zPrimeA(:,1:2)=zPrime;
        zPrimeB=zPrimeA;
        %%%rescale z values
        %Integral A, dt1
        zPrimeAL=tauPrime/2;
        zPrimeAU=min(1/gammaConst/2,tauPrime);
        zPrimeWidthA=(zPrimeAU-zPrimeAL);
        zPrimeA(:,1)=zPrimeWidthA.*zPrimeA(:,1)+zPrimeAL;
        zPrimeVolA=zPrimeVolA.*zPrimeWidthA;
        %Integral A, dt2
        zPrimeAL=tauPrime-zPrimeA(:,1);
        zPrimeAU=zPrimeA(:,1);
        zPrimeWidthA=(zPrimeAU-zPrimeAL);
        zPrimeA(:,2)=zPrimeWidthA.*zPrimeA(:,2)+zPrimeAL;
        zPrimeVolA=zPrimeVolA.*zPrimeWidthA;
        %volume part
        gammaIndA=(zPrimeA(:,1)+(i+1)*zPrimeA(:,2))<(1/gammaConst);
        epsilonIndA=zPrimeA(:,2)>epsilonPrime;
        icIndA=(zPrimeA(:,1)+gammaConst*tauPrime*zPrimeA(:,2))>tauPrime;
        zPrimeVolA=zPrimeVolA.*epsilonIndA.*gammaIndA;
        zPrimeVol_ICA=zPrimeVolA.*icIndA;
        zPrimeCumA=zPrimeA(:,1)+zPrimeA(:,2);
        
        if tauPrime>(1/gammaConst/(2+0))
            %Integral B, dt1
            zPrimeBL=1/gammaConst/2;
            zPrimeBU=tauPrime;
            zPrimeWidthB=(zPrimeBU-zPrimeBL);
            zPrimeB(:,1)=zPrimeWidthB.*zPrimeB(:,1)+zPrimeBL;
            zPrimeVolB=zPrimeVolB.*zPrimeWidthB;
            %Integral B, dt2
            zPrimeBL=tauPrime-zPrimeB(:,1);
            zPrimeBU=(1/gammaConst-zPrimeB(:,1))/(1+i);
            zPrimeWidthB=(zPrimeBU-zPrimeBL);
            zPrimeB(:,2)=zPrimeWidthB.*zPrimeB(:,2)+zPrimeBL;
            zPrimeVolB=zPrimeVolB.*zPrimeWidthB;
            %volume part
            gammaIndB=(zPrimeB(:,1)+(i+1)*zPrimeB(:,2))<(1/gammaConst);
            epsilonIndB=zPrimeB(:,2)>epsilonPrime;
            zPrimeVolB=zPrimeVolB.*epsilonIndB.*gammaIndB;
            icIndB=(zPrimeB(:,1)+tauPrime*zPrimeB(:,2))>tauPrime;
            zPrimeVol_ICB=zPrimeVolB.*icIndB;
            zPrimeB(:,(nComb+1):end)=repmat(zPrimeB(:,2),1,i);
            
        else
            zPrimeVolB=0; zPrimeVol_ICB=0;
        end
        zPrimeA(:,(nComb+1):end)=repmat(zPrimeA(:,2),1,i);
        
        %create/sample dQt1dt2 function
        %using derivative function
        dQt1dt2A=funQndt1dt2(zPrimeA,betaConst,gammaConst,etaValues,(i+2));
        if tauPrime>(1/gammaConst/(2+0))
            dQt1dt2B=funQndt1dt2(zPrimeB,betaConst,gammaConst,etaValues,(i+2));
        else
            dQt1dt2B=0;
        end
        
        %generate integral kernel and estimate integral
        kernelInt=RviRand.*(dQt1dt2A.*zPrimeVolA+dQt1dt2B.*zPrimeVolB)./funIS; %integral kernerl
        mu_hat_n=mean(kernelInt); % perform (q)MC step and add it to the previous result
        In=funIn(betaConst,nMoment,x); %calculate In integral eq. (12) in [1]
        mu_n=In*mu_hat_n/nMoment*factorial(nMoment); %density of Mn
        ithTerm=(-1)^i*mu_n/factorial(i);
        probDeltaMCComb(k,(i+1):end)=probDeltaMCComb(k,(i+1):end)+ithTerm;
        
        %%_IC section
        %generate integral kernel and estimate integral
        kernelInt=RviRand.*(dQt1dt2A.*zPrimeVol_ICA+dQt1dt2B.*zPrimeVol_ICB)./funIS; %integral kernerl
        mu_hat_n=mean(kernelInt); % perform (q)MC step and add it to the previous result
        In=funIn(betaConst,nMoment,x); %calculate In integral eq. (12) in [1]
        mu_n=In*mu_hat_n/nMoment*factorial(nMoment); %density of Mn
        ithTerm=(-1)^i*mu_n/factorial(i);
        probDeltaMC_IC(k,(i+1):end)=probDeltaMC_IC(k,(i+1):end)+ithTerm;
        
    end
end

%plotting section
tauValues=tauPrimeValues./(1-gammaConst*tauPrimeValues);
tauValuesdB=10*log10(tauValues);

figure;
set(gcf,'DefaultLineLineWidth',2);
set(gca, 'FontSize', 14);
plot(tauValuesdB,probDeltaMCComb(:,end),tauValuesdB,probDeltaMC_IC(:,end),'--');grid;
title('Integration Method');
xlabel('$\tau (dB)$','Interpreter','latex');
ylabel('$\Delta(\tau,\epsilon)$','Interpreter','latex');
legend('Signal Combination','Interference Cancellation');
title('');
