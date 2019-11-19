% Calculates k-coverage probability for a single-tier cellular network
% with log-normal shadowing based on a Poisson model outlined in [1]
% For SIR results (ie noise N=0) results are scale-invariant (independent
% of lambda)
% A variation of this code was used to produce Figures 1 and 2 in [2]
% Author: H.P. Keeler, Inria Paris/ENS, 2014
%
% References:
% [1] H.P. Keeler, B. Błaszczyszyn and M. Karray,
% 'SINR-based coverage probability in cellular networks with arbitrary 
% shadowing', ISIT, 2013
% [2] B. Błaszczyszyn  and H.P. Keeler, 'Studying the SINR process in
% Poisson networks by using its factorial moment measures', submitted,
% 2013
% [3] Błaszczyszyn  and H.P. Keeler, 'Equivalence and comparison of random
% heterogeneous networks', presented at WDN-CN, 2013
% [4] Dhillon, Ganti, Adnrews, Baccelli, 'Modeling and analysis of K-tier
% downlink heterogeneous cellular networks', JSAC, 2012
%
% NOTE: If you use this code in published work, please cite paper[2] by
% Błaszczyszyn and Keeler, as listed above.



clear all; close all; clc;

lambda=1; %base station density

%Power terms
P=1; 

betaConst=3; %path-loss exponent  
K=P.^(-betaConst);

%log normal parameters
sigmDb=10;
sigma=sigmDb/10*log(10);
ESTwoBeta=exp(sigma^2*(2-betaConst)/betaConst^2);

%model constant - incorporates model parameters
a=lambda*pi*ESTwoBeta/K^2; %equation (6) in [1]
x=W*a^(-betaConst/2); %this constant represents all the network parameters except the path-loss exponent

%noise parameters
N=0;10^(-96/10)/1000;
P=10^(62.2/10)/1000;
W=N/P;

%SINR threshold values 
tMinDb=-10;tMaxDb=4;
tValuesDb=(tMinDb:tMaxDb)'; %values in dB
tValues=10.^(tValuesDb/10);
tNumb=length(tValues);

k=1; %coverage number
%analytic/integration section
numbMC=10^3;
PCov1=funProbCov(tValues,betaConst,x,numbMC,k);
k=2; %coverage number
PCov2=funProbCov(tValues,betaConst,x,numbMC,k);
k=3; %coverage number
PCov3=funProbCov(tValues,betaConst,x,numbMC,k);


%plotting section

Pn1=1-PCov1;
Pn2=1-PCov2;
Pn3=1-PCov3;

%create suitable label
if W==0
    legendLabel='SIR';
else 
    legendLabel='SINR';
end

figure;
set(gcf,'DefaultLineLineWidth',2);
set(gca, 'FontSize', 14);
plot(tValuesDb,PCov1,'-o',tValuesDb,PCov2,tValuesDb,PCov3,'-x');grid;
legend('k=1','k=2','k=3','Location','NorthEast');
xlabel('$\tau (dB)$','Interpreter','latex'); 
ylabel('$\mathcal{P}^{(k)}(\tau)$','Interpreter','latex');
axis([tValuesDb(1),tValuesDb(end),0,1]);

