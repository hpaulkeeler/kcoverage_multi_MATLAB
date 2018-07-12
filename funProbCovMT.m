% funProbCovMT calculates and returns SINR-based k-th coverage probability
% under arbitrary shadowing in a multi-tier cellular network based on a
% model outlined in [3,4] with calculation methods developed in [3]
% PCovMT=funProbCovMT(tValues,k,betaConst,KValues,lambdaValues,
% ... ESTwoBetaValues,W,gammaConst,numbMC)
% PCovMT is the k-coverage probability
% tValues are the SINR threshold values. tValues can be a vector
% k is coverage number ie maximum number of connected base stations
% betaConst is the path-loss exponent for the entire network
% KValues is the path-loss constant for each tier
% lambdaValues is the base station density for each tier
% ESTwoBetaValues are the E(S^(2/beta)) fading/shadowing moments for tiers
% Length of tValues is the number of tiers
% Lengths of tValues, KValues, lambdaValues and ESTwoBetaValues must match
% numbMC is number of quasi-random points for the numerical integratopm
% betaConst, W, gamma  and numbMC are scalars.
%
% Author: H.P. Keeler, Inria/ENS Paris, 2018
%
% References:
% [1] H.P. Keeler, B. Błaszczyszyn and M. Karray,
% 'SINR-based coverage probability in cellular networks with arbitrary 
% shadowing', ISIT, 2013
% [2] B. Błaszczyszyn  and H.P. Keeler, 'Studying the SINR process in
% Poisson networks by using its factorial moment measures' , IEEE TOIT,
% 2015
% [3] Błaszczyszyn  and H.P. Keeler, 'Equivalence and comparison of random
% heterogeneous networks', presented at WDN-CN, 2013
% [4] Dhillon, Ganti, Adnrews, Baccelli, 'Modeling and analysis of K-tier
% downlink heterogeneous cellular networks', JSAC, 2012


function PCovMT=funProbCovMT(tValues,k,betaConst,KValues,lambdaValues,ESTwoBetaValues,W,gammaConst,numbMC)

m=length(tValues); %number of tiers
%check vector lengths match
m1=length(KValues); m2=length(lambdaValues); m3=length(ESTwoBetaValues);
check_length=sum(abs([m1,m2,m3]-m)); %equal to zero for equal lengths
if check_length>0
    error('The lengths of tValues, KValues, lambdaValues and ESTwoBetaValues do not match!');
else
    %maxium number of symmetric terms
    nValues=max(1,ceil(1./(gammaConst*tValues)));
    nMax=max(nValues);
    %Multi-tier section
    lambdaValues_star=lambdaValues./KValues.^2.*ESTwoBetaValues;
    lambda_star=sum(lambdaValues_star); %single-tier base station density
    a=lambda_star*pi;%Uses equivalent network: S=1 so E(S^(2/beta)=1), K=1
    x=(W/gammaConst)*a.^(-betaConst/2);
    PCovMT=0;
    SnMT=zeros(1,nMax);
    for n=k:nMax
        SnMT(n)=0;
        permNumb=m^n; %number of permutations
        tValuesPerm=zeros(permNumb,n); lambdaValuesPerm=zeros(permNumb,n);
        for i=1:permNumb
            %create permutations of threshold and lambda values
            for j = n:-1:1
                tmp = repmat(tValues',m^(n-j),m^(j-1));
                tValuesPerm(:,j) = tmp(:);
                tmp = repmat(lambdaValues_star',m^(n-j),m^(j-1));
                lambdaValuesPerm(:,j) = tmp(:);
            end
            t_temp=tValuesPerm(i,:);
            t_temp_prime=t_temp./(1+gammaConst*t_temp);
            if (gammaConst*sum(t_temp_prime))<1
                lambda_temp=lambdaValuesPerm(i,:);
                t_hat=gammaConst*t_temp_prime./(1-gammaConst*sum(t_temp_prime));
                In=funIn(betaConst,n,x); %eq. (12) in [1]
                %calcualtes Jn integral for the multi-tier setting
                JnMT=funJnMT(t_hat,betaConst,length(t_hat),numbMC,1);
                Mn=prod(t_hat.^(-2/betaConst)).*In.*JnMT;
                SN_perTierPerm=Mn.*prod(lambda_temp/lambda_star);
            else
                SN_perTierPerm=0;
            end
            SnMT(n)=SnMT(n)+SN_perTierPerm;
        end
        SnMT(n)=SnMT(n);
        PCovMT=(-1)^(n-k)*nchoosek(n-1,k-1)*SnMT(n)+PCovMT; %eq. (8) in [1]
    end
end
