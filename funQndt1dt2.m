% Calculates the second-order "partial density" of the n th factorial
% moment measure of the STINR processas outlined Section C of the
% Appendix in [2]. This is used for obtaining the joint distribution
% of the first and second order statistics of the STINR process.  
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


 function Qn_12=funQndt1dt2(tPrimeValues,betaConst,gammaConst,etaValues,n)
  
 tPrime=tPrimeValues;
 eta=etaValues;


%some constants to ease calculation
 alphaConst=2/betaConst;
b=n*(alphaConst+1);

%%% A function section
A=(gammaConst)^(-b)*(1-gammaConst*sum(tPrime,2)).^(b); %A function definition
%derivatives of A function
A_1=-b*(gammaConst)^(-b+1)*(1-gammaConst*sum(tPrime,2)).^(b-1);
A_2=A_1;
A_12=b*(b-1)*(gammaConst)^(-b+2)*(1-gammaConst*sum(tPrime,2)).^(b-2);

tPrimeAlpha=tPrime.^(alphaConst);
tPrimeAlpha1=tPrime.^(alphaConst+1);

%%% B function section
B=1./prod(tPrimeAlpha1+bsxfun(@times,eta./gammaConst.*tPrimeAlpha,...
    (1-gammaConst.*sum(tPrime,2))),2); %B function definition
cumB=1./cumprod(tPrimeAlpha1+bsxfun(@times,eta./gammaConst.*tPrimeAlpha,...
    (1-gammaConst*sum(tPrime,2))),2);


%sum remaing t values ie not t1 and t2, which willb e used in P function
%definition
tPrimeSum=sum(tPrime(:,3:end),2);

%D function definitions
D1=tPrimeAlpha1(:,1)...
    +eta(:,1)./gammaConst.*tPrimeAlpha(:,1).*(...
    1-gammaConst.*tPrimeSum);
D2=tPrimeAlpha1(:,2)...
    +eta(:,2)./gammaConst.*tPrimeAlpha(:,2).*(...
    1-gammaConst.*tPrimeSum);

%derivatives of D function
D1_1=(alphaConst+1)*tPrimeAlpha(:,1)...
    +eta(:,1)/gammaConst.*alphaConst.*tPrime(:,1).^(alphaConst-1).*(...
    1-gammaConst.*tPrimeSum);
D2_2=(alphaConst+1)*tPrimeAlpha(:,2)...
     +eta(:,2)/gammaConst.*alphaConst.*tPrime(:,2).^(alphaConst-1).*(...
    1-gammaConst.*tPrimeSum);

%B function definitions
B1=(D1-eta(:,1).*(tPrimeAlpha1(:,1)+tPrimeAlpha(:,1).*tPrime(:,2))).^(-1);
B2=(D2-eta(:,2).*(tPrimeAlpha1(:,2)+tPrimeAlpha(:,2).*tPrime(:,1))).^(-1);

%derivatives of B functions
B1_1=-1*(D1_1-eta(:,1).*((alphaConst+1).*tPrimeAlpha(:,1)+alphaConst*tPrime(:,1).^(alphaConst-1).*tPrime(:,2))).*B1.^2;
B2_2=-1*(D2_2-eta(:,2).*((alphaConst+1).*tPrimeAlpha(:,2)+alphaConst*tPrime(:,2).^(alphaConst-1).*tPrime(:,1))).*B2.^2;

B1_12=eta(:,1).*alphaConst.*tPrime(:,1).^(alphaConst-1).*(B1).^2-...
    2*(D1_1-eta(:,1).*((alphaConst+1).*tPrimeAlpha(:,1)+alphaConst*tPrime(:,1).^(alphaConst-1).*tPrime(:,2)))...
    .*(eta(:,1).*tPrime(:,1).^alphaConst).*B1.^3;
B2_12=eta(:,2).*alphaConst.*tPrime(:,2).^(alphaConst-1).*(B2).^2-...
    2*(D2_2-eta(:,2).*((alphaConst+1).*tPrimeAlpha(:,2)+alphaConst*tPrime(:,2).^(alphaConst-1).*tPrime(:,1)))...
    .*(eta(:,2).*tPrime(:,2).^alphaConst).*B2.^3;

B2_1=eta(:,2).*tPrimeAlpha(:,2).*B2.^2;
B1_2=eta(:,1).*tPrimeAlpha(:,1).*B1.^2;

if n>=3
    %%% section for calculate product function P    
    P=B./B1./B2; %P function definition
    %cumP=cumB(:,2:end)./B1./B2; %acculative version of P function
    cumP=bsxfun(@rdivide,cumB(:,2:end),B1.*B2); %acculative version of P function
    
    
    %calculate B functions
    Bi=cumP(:,2:end)./cumP(:,1:end-1);    
    sumBi=sum(bsxfun(@times,Bi,eta(:,3:end).*tPrimeAlpha(:,3:end)),2);    
    %derivatives of P function
    P_1=P.*sumBi;
    P_2=P_1;
    P_12=P_2.*(sumBi)+P.*sum((bsxfun(@times,Bi,(eta(:,3:end).*tPrimeAlpha(:,3:end))).^2),2);
    B_12=B1_12.*B2.*P+B1_1.*B2_2.*P+B1_2.*B2_1.*P...
        +B1.*B2_12.*P+B1.*B2_1.*P_2+B1.*B2_2.*P_1...
        +B1_1.*B2.*P_2+B1_2.*B2.*P_1+B1.*B2.*P_12;
    B_1=B1_1.*B2.*P+B1.*B2_1.*P+B1.*B2.*P_1;
    B_2=B1_2.*B2.*P+B1.*B2_2.*P+B1.*B2.*P_2;
else
    %since P=1 for n=2, derivatives simplify
    B_12=B1_12.*B2+B1_1.*B2_2+B1_2.*B2_1+B1.*B2_12;    
    B_1=B1_1.*B2+B1.*B2_1;
    B_2=B1_2.*B2+B1.*B2_2;
end

%%% hj function section
hj=gammaConst*bsxfun(@rdivide,tPrime,(1-gammaConst*sum(tPrime,2))); %hj function definition
%derivatives of h function
h1_1=gammaConst*bsxfun(@rdivide,(1-gammaConst*(sum(tPrime,2)-tPrime(:,1))),(1-gammaConst*sum(tPrime,2)).^2);
h2_2=gammaConst*bsxfun(@rdivide,(1-gammaConst*(sum(tPrime,2)-tPrime(:,2))),(1-gammaConst*sum(tPrime,2)).^2);

hj_1=gammaConst^2*bsxfun(@rdivide,tPrime,(1-gammaConst*sum(tPrime,2)).^2);  %j > 2
hj_2=gammaConst^2*bsxfun(@rdivide,tPrime,(1-gammaConst*sum(tPrime,2)).^2); %%j > 2
hj_12=2*gammaConst^3*bsxfun(@rdivide,tPrime,(1-gammaConst*sum(tPrime,2)).^3); %j > 2

h1_12=gammaConst^2.*(-1+2*(1-gammaConst*(sum(tPrime,2)-tPrime(:,1)))...
    .*(1-gammaConst*sum(tPrime,2)).^(-1))./(1-gammaConst*sum(tPrime,2)).^2;
h2_12=gammaConst^2.*(-1+2*(1-gammaConst*(sum(tPrime,2)-tPrime(:,2)))...
    .*(1-gammaConst*sum(tPrime,2)).^(-1))./(1-gammaConst*sum(tPrime,2)).^2;

%%% calculate derivative of Q0 using A and B functions
Q0=A.*B;
Q0_1=A_1.*B+A.*B_1; 
Q0_2=A_2.*B+A.*B_2;
Q0_12=A_12.*B+A_1.*B_2+A_2.*B_1+A.*B_12;

%%%calculate sum_j(hj QO) using h_i and Q derivatives
sumhj=sum(hj,2);
sumhj_2=sum(hj_2,2)-hj_2(:,2)+h2_2;
sumhj_1=sum(hj_1,2)-hj_2(:,1)+h1_1;
sumhj_12=sum(hj_12,2)-hj_12(:,1)-hj_12(:,2)+h1_12+h2_12;
sum_Qi_12=Q0_12.*sumhj+Q0_1.*sumhj_2+Q0_2.*sumhj_1+Q0.*sumhj_12;


%%% calculate sum_i(hj QO) using h_i and Q derivatives
Qn_12=Q0_12+sum_Qi_12;
