% Calculates and plots SINR-based k-coverage probability in a multi-tier
% network and a single-tier equivalent network. The single-tier equivalent 
% has the same propation processes with random SINR threshold values
% averaged; see [3] for more details on equivalent Poisson networks.
% Integration method based on Sobol points - outlined in [2].
% A variation of this code was used to produce Figures 3 to 5 in [2]
% Warning: Be careful when setting tMinDb1/tMinDb2 as integration method 
% may take a long time for values less than -10 dB.
% Increase the variable numbMC, number of sample points, for more accuracy 
% in the integration method but usually 2^10 is sufficient.
% 
% Author: H.P. Keeler, Inria Paris/ENS, 2014
%
% References
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

% Parameters for Fig. 7 in [4]
% P2=1; P1=100*P2; power terms, N=0; no noise term (so scale-invariant)
% lambda2= 2*lambda1; densities
% betaConst=3 %path-loss exponent
% tau2=1 dB; other

clear all; clc; close all;

k=2; %coverage number number

betaConst=3; %path-loss exponent

%Power terms
P2=1; P1=100*P2; PValues=[P1;P2];

%base station density parameters
lambda1=1; lambda2=2*lambda1;
lambdaValues=[lambda1,lambda2]; lambdaValues=lambdaValues(:);
onesVector=ones(size(lambdaValues));

m=length(lambdaValues); %number of tiers

%path-loss constant
KValues=PValues.^(-betaConst); KValues=KValues(:);
%single-tier values
PST=1; %single-tier equivalent power term
KST=1;%single-tier equivalent path-loss constant 
lambdaST=sum(lambdaValues.*PValues.^(2/betaConst)); %single-tier base station density;

%inteference cancellation factor
gammaConst=1;

%exponential variables
lambdaE=1;
ESTwoBetaValues=lambdaE^(2/betaConst)*gamma(2/betaConst+1)*onesVector;

%noise parameters
N=0;%10^(-96/10)/1000;
P=10^(62.2/10)/1000;
W=N/P;

%SINR threshold values
tMinDb1=-10;tMaxDb1=10; %Dhillon et al. [3] range
tValuesDb1=(tMinDb1:tMaxDb1)'; %values in dB
tValues1=10.^(tValuesDb1/10);
tNumb1=length(tValues1);

tValuesDb2=1; %values in dB
tValues2=10.^(tValuesDb2/10);

%analytic/integration section
numbMC=2^8; %number of sample points
%find equivalent lambda values
lambdaValues_star=lambdaValues./KValues.^2.*ESTwoBetaValues;
lambda_star=sum(lambdaValues_star); %single-tier base station density
a=lambda_star*pi;%S=1 so E(S^(2/beta)=1), K=1
x=(W/gammaConst)*a.^(-betaConst/2);

%Single-tier
ESTwoBetaST=sum(lambdaValues.*ESTwoBetaValues)/lambda_star; %spatial average
lambda_starST=lambdaST./KST.^2.*ESTwoBetaST; %single-tier base station density
aST=lambda_starST*pi;%S=1 so E(S^(2/beta)=1), K=1
xST=(W/gammaConst)*aST.^(-betaConst/2);

%Simulation section
simNumb=10^3; %number of simulations
diskRadius=10; %radius of simulation disk region (has to be larger when fading is incorporated)
PCov=zeros(tNumb1,1); PCovST=zeros(tNumb1,1);
%vary t_1 (first tier) values, and hold t_2 (second-tier) values fixed.
for i=1:tNumb1
    tValues=[tValues1(i);tValues2];
    PCov(i)=funProbCovMT(tValues,k,betaConst,KValues,lambdaValues,ESTwoBetaValues,W,gammaConst,numbMC);
    tST=sum(tValues.*lambdaValues.*PValues.^(2/betaConst))/lambdaST;
    PCovST(i)=funProbCov(tST,betaConst,W*aST^(-betaConst/2),numbMC,k);
end

%plotting section
Pn=1-PCov;
PnSim=1-PCovST;

%create suitable label
if W==0
    legendLabel='SIR';
else
    legendLabel='SINR';
end

figure;
set(gcf,'DefaultLineLineWidth',2);
set(gca, 'FontSize', 14);
plot(tValuesDb1,PCov,'-o',tValuesDb1,PCovST,'-');grid;
legend('Two-tier','Single-tier','Location','NorthEast')
xlabel('$\tau_1 (dB)$','Interpreter','latex');
ylabel('$\mathcal{P}^{(k)}(\tau_1,\tau_2)$','Interpreter','latex');
axisValues=[min(tValuesDb1), max(tValuesDb1), 0, 1];axis(axisValues);
title('k-coverage Probability');
