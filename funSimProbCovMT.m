% funSimSnMT estimates and returns the calculates and returns k-th coverage
% probability under arbtirary shadowing based on a simulation (on a disk)
% of the model in [2], which extends models in [1] and [4]
% simPCovMT=funSimProbCovMT(tValues,k,betaConst,KValues,lambdaValues,
% ESTwoBetaValues,W,gammaConst,simNumb,diskRadius)
% simPCovMT is the n-th symmetric sum
% tValues are the SINR threshold values. tValues can be a vector
% k is the coverage ie maximum number of connected base stations
% betaConst is the path-loss exponent for the entire network
% KValues is the path-loss constant for each tier
% lambdaValues is the base station density for each tier
% ESTwoBetaValues are the E(S^(2/beta)) fading/shadowing moments the tiers
% Length of tValues is the number of tiers
% Lengths of tValues, KValues, lambdaValues and ESTwoBetaValues  must match
% simPCovMT is the k-coverage probability
% diskRadius is radius of network disk
% betaConst, W, gamma  and numbMC are scalars.
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
%
% NOTE: If you use this code in published work, please cite paper[2] by
% Błaszczyszyn and Keeler, as listed above.

                     
function simPCovMT=funSimProbCovMT(tValues,k,betaConst,KValues,lambdaValues,ESTwoBetaValues,W,gammaConst,simNumb,diskRadius)

m=length(tValues); %number of tiers of base stations (BSs)
%check vector lengths match
m1=length(KValues); m2=length(lambdaValues); m3=length(ESTwoBetaValues);
check_length=sum(abs([m1,m2,m3]-m)); %equal to zero for equal lengths
if check_length>0
    error('The lengths of tValues, KValues, lambdaValues and ESTwoBetaValues do not match!');
else
    diskArea=pi*diskRadius^2;  %area of simulation region
    coveredNumbk=0;
    %rescale lambda - see foot note 5 in [1]
    lambdaSimValues=lambdaValues.*ESTwoBetaValues; 
    for i=1:simNumb
        
        randNumb=poissrnd(lambdaSimValues*diskArea); %# of BSs per tier
        randNumbTotal=sum(randNumb); %total # of BSs in network
        tPerBS=zeros(randNumbTotal,1);
        KPerBS=zeros(randNumbTotal,1);
        firstIndex=1; lastIndex=randNumb(1);
        for l=1:length(lambdaValues)
            %done this way to preallocate vectors
            tPerBS(firstIndex:lastIndex)=tValues(l);
            KPerBS(firstIndex:lastIndex)=KValues(l);
            if l<m
                firstIndex=lastIndex+1;
                lastIndex=firstIndex+randNumb(l+1)-1;
            end
        end
        %random distances from the typical node
        rRand=diskRadius*sqrt(rand(randNumbTotal,1)); %uniform in cartesion, not polar coordinates
        %shadowing distribution can be constant if lambda is rescaled - see [1]
        shadowRand=ones(randNumbTotal,1);
        signalRand=shadowRand.*(KPerBS.*rRand).^(-betaConst);
        interferTotal=sum(signalRand); %total inteference in network
        SINR=signalRand./(gammaConst*(interferTotal-signalRand)+W); %calculate SINR for each node in the network
        
        %counts how many nodes are exactly k or more connected/covered
        if sum(SINR>tPerBS)>=k
            coveredNumbk=coveredNumbk+1;
        end        
    end
    
    simPCovMT=coveredNumbk/simNumb;   %calculate k-coverage probability
end
