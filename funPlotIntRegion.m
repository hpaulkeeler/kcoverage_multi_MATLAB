% Plots the integration region for two signal combination
% model ie illustrates the region over which to integrate 
% the joint density of the orders statistics of the STINR process  
% This code is not actually used to calculate anything, but useful 
% for visualizing and understanding the integral region

% Author: H.P. Keeler, Inria/ENS, Paris and University of Melbourne, 2018
%
% References
% [1] H.P. Keeler, B. Błaszczyszyn and M. Karray,
% 'SINR-based coverage probability in cellular networks with arbitrary
% shadowing', ISIT, 2013
% [2] B. Błaszczyszyn  and H.P. Keeler, 'Studying the SINR process in
% Poisson networks by using its factorial moment measures', IEEE TOIT,
% 2015

function funPlotIntRegion(gammaConst,epsilon,tauPrime,i)

%close all; clear all;clc;
% gammaConst=1;
% tauPrime=.6;
% epsilon=0.00001;

tNumb=20;
t1=linspace(0, 1/gammaConst,tNumb);

%first line
%t1>t2
line1=t1;

%second line
%t1+(t2+i)<1/gammaConst
line2=(1/gammaConst-t1)/(i+1);

%third line
%t1+t2>tauPrime
t1_3=linspace(0, tauPrime,tNumb);
line3=tauPrime-t1_3;

%fourth line
%t2>epsilon
line4=epsilon*ones(1,tNumb);

figure;
plot(t1,line1,t1,line2,t1_3,line3,t1,line4);
legend('Line 1','Line 2','Line 3','Line 4','Location','North');
axis equal; axis tight;
hold on;
%plot tauPrime square
plot(tauPrime*ones(1,tNumb),linspace(0,tauPrime,tNumb)); 
plot(linspace(0,tauPrime,tNumb),tauPrime*ones(1,tNumb));

%plot interection points
scatter(tauPrime/2,tauPrime/2,'o'); %t1=t2=tauPrime/2
scatter(1/gammaConst/(2+i),1/gammaConst/(2+i),'*'); %t1=t2=1/gammaConst/2

titleStr1=sprintf('i = %d',i);
titleStr2=sprintf('tauPrime = %0.5g',tauPrime);
titleStr3=sprintf('epsilon = %0.5g',epsilon);
titleStr=[titleStr1, ', ',titleStr2, ', ',titleStr3];
title(titleStr);

end
