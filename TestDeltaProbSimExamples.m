% Estimates delta probability for 2-base station signal combination
% and second-strongest signal removed under arbitrary shadowing/fading
% in a multi-tier cellular network based on a Poisson
% model outlined in [2] (which is similar to [4]).
% For SIR results (ie noise N=0) results are scale-invariant (independent
% of base station desnity lambda)
% In the simulation, the base stations are poistioned on a disk
% A variation of this code was used to produce Figures 6 to 11 in [2]
%
% Author: H.P. Keeler, Inria Paris/ENS, 2014
%
% References:
% [1] H.P. Keeler, B. Błaszczyszyn and M. Karray,
% 'SINR-based coverage probability in cellular networks with arbitrary
% shadowing', ISIT, 2013
% [2] B. Błaszczyszyn  and H.P. Keeler, 'Studying the SINR process in
% Poisson networks by using its factorial moment measures' , submitted,
% 2013
% [3] Błaszczyszyn  and H.P. Keeler, 'Equivalence and comparison of random
% heterogeneous networks', presented at WDN-CN, 2013
% [4] Dhillon, Ganti, Adnrews, Baccelli, 'Modeling and analysis of K-tier
% downlink heterogeneous cellular networks', JSAC, 2012
clear all; close all; clc;

simNumb=10^4; %% number of network simulations

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
epsilonPrime=0.1; %lower threshold value
tauNumb=25;
tauPrimeValues=linspace(epsilonPrime,.99,tauNumb);

%noise parameters
N=0;10^(-96/10)/1000;
P=10^(62.2/10)/1000;
W=N/P;

%%% Simulation Section %%%
%(uniformly) randomly places nodes on a disk of radius diskRadius
%simulation parameters
diskRadius=50;
diskArea=pi*diskRadius^2;
numbDeltaComb=zeros(size(tauPrimeValues));
numbDeltaic=zeros(size(tauPrimeValues));

%rescale lambda - see foot note 5 in [1]
lambdaSim=lambda.*ESTwoBeta;
for i=1:simNumb
    randNumb=poissrnd(lambdaSim*diskArea);
    randNumbTotal=sum(randNumb);
    %random distances from the typical node
    rRand=diskRadius*sqrt(rand(randNumbTotal,1)); %uniform in cartesion, not polar coordinates
    %shadowing distribution can be constant if lambda is rescaled - see [1]
    shadowRand=ones(randNumbTotal,1);
    signalRand=shadowRand.*(K.*rRand).^(-betaConst);
    interferTotal=sum(signalRand); %total inteference in network
    SINR=signalRand./(gammaConst*(interferTotal)+W); %calculate modified SINR for each node in the network
    %counts how many nodes are exactly k or more connected/covered
    sortedSINR=sort(SINR,'descend');
    coordSINR=sum(sortedSINR(1:nComb));
    for k=1:tauNumb
        tauPrime=tauPrimeValues(k);
        if (coordSINR>tauPrime)&&(sortedSINR(2)>epsilonPrime)&&(sortedSINR(1)<tauPrime)
            numbDeltaComb(k)=numbDeltaComb(k)+1;
        end
        if (sortedSINR(1)+gammaConst*tauPrime*sortedSINR(2)>tauPrime)&&(sortedSINR(2)>epsilonPrime)&&(sortedSINR(1)<tauPrime)
            numbDeltaic(k)=numbDeltaic(k)+1;
        end
    end
end
probDeltaSimComb=numbDeltaComb/simNumb;
probDeltaSimic=numbDeltaic/simNumb;

%plotting section
tauValues=tauPrimeValues./(1-gammaConst*tauPrimeValues);
tauValuesdB=10*log10(tauValues);

figure;
set(gcf,'DefaultLineLineWidth',2);
set(gca, 'FontSize', 14);
plot(tauValuesdB,probDeltaSimComb,'*',tauValuesdB,probDeltaSimic,'o');grid;
title('Simulation Method');
xlabel('$\tau (dB)$','Interpreter','latex');
ylabel('$\Delta(\tau,\epsilon)$','Interpreter','latex');
legend('Signal Combination','Interference Cancellation');
shg;
axis tight;

