% Calculates and plots SINR-based k-coverage probability in a multi-tier
% network by both simulating a network model and using an integration method 
% outlined [2]
%
% Increase the variable numbMC, number of sample points, for more accuracy 
% in the integration method but usually 2^10 is sufficient.
%
% Author: H.P. Keeler, Inria/ENS, Paris and University of Melboure, 2018
%
% References
% [1] H.P. Keeler, B. Błaszczyszyn and M. Karray,
% 'SINR-based coverage probability in cellular networks with arbitrary
% shadowing', ISIT, 2013
% [2] B. Błaszczyszyn  and H.P. Keeler, 'Studying the SINR process in
% Poisson networks by using its factorial moment measures', IEEE TOIT,
% 2015


clear all; close all; clc;


numbMC=2^8; %number of sample points for integral estimate

k=1 ; %coverage number number

%base station density parameters
lambdaValues=[3,2.5]; lambdaValues=lambdaValues(:);
lambda=sum(lambdaValues); %single-tier base station density
onesVector=ones(size(lambdaValues));

m=length(lambdaValues); %number of tiers

%path-loss parameters
betaConst=3.8; %path-loss exponent
K=1; %path-loss constant
KValues=K*[2,17]; KValues=KValues(:); %vector of path-loss values

%inteference cancellation factor
gammaConst=1;

%log normal parameters
sigmDb=10;
sigma=sigmDb/10*log(10);
ESTwoBetaValues=exp(sigma^2*(2-betaConst)/betaConst^2)*onesVector;

%noise parameters
N=10^(-96/10)/1000;
P=10^(62.2/10)/1000;
W=0; %W=N/P;


%SINR threshold values
tMinDb1=-2;tMaxDb1=5;
tValuesDb1=(tMinDb1:tMaxDb1)'; %values in dB
tValues1=10.^(tValuesDb1/10);
tNumb1=length(tValues1);

tValues2=1;

%find equivalent lambda values
lambdaValues_star=lambdaValues./KValues.^2.*ESTwoBetaValues;
lambda_star=sum(lambdaValues_star); %single-tier base station density
a=lambda_star*pi;%S=1 so E(S^(2/beta)=1), K=1
x=(W/gammaConst)*a.^(-betaConst/2); 

PCov=zeros(1,tNumb1); PCovSim=zeros(1,tNumb1);
%vary t_1 (first tier) values, and hold t_2 (second-tier) values fixed.

for i=1:tNumb1
    tValues=[tValues1(i),tValues2];
    PCov(i)=funProbCovMT(tValues,k,betaConst,KValues,lambdaValues,ESTwoBetaValues,W,gammaConst,numbMC);
    
    %Simulation section
    simNumb=10^3; %number of simulations
    diskRadius=10; %radius of simulation disk region (has to be larger when fading is incorporated)
    PCovSim(i)=funSimProbCovMT(tValues,k,betaConst,KValues,lambdaValues,ESTwoBetaValues,W,gammaConst,simNumb,diskRadius);
    
end

%plotting section
Pn=1-PCov;
PnSim=1-PCovSim;

%create suitable label
if W==0
    legendLabel='SIR';
else
    legendLabel='SINR';
end

figure;
set(gcf,'DefaultLineLineWidth',2);
set(gca, 'FontSize', 14);
plot(tValuesDb1,PCovSim,'o',tValuesDb1,PCov);grid;
legend(['Simulation ', legendLabel],['Integration ', legendLabel],'Location','NorthWest')
xlabel('$\tau_1 (dB)$','Interpreter','latex');
ylabel('$\mathcal{P}^{(k)}(\tau_1,\tau_2)$','Interpreter','latex');
