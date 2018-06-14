% Logan Prust - Aer E 311 - Homework 2 - Problem 5

clear,clc;

gamma=1.4;

% get input parameters
M1=input('Upstream Mach Number: ');
P1=input('Upstream Pressure: ');
rho1=input('Upstream Density: ');
T1=input('Upstream Temperature: ');
Pt1=input('Upstream Total Pressure: ');

% calculate output parameters
M2=sqrt((M1^2+2/(gamma-1))/(2*gamma/(gamma-1)*M1^2-1));
P2=P1*(1+gamma*M1^2)/(1+gamma*M2^2);
T2=T1*(1+(gamma-1)/2*M1^2)/(1+(gamma-1)/2*M2^2);
rho2=rho1*M1/M2*sqrt((1+(gamma-1)/2*M2^2)/(1+(gamma-1)/2*M1^2));
Pt2=Pt1*(1+gamma*M1^2)/(1+gamma*M2^2)*((1+(gamma-1)/2*M2^2)/(1+(gamma-1)/2*M1^2))^(gamma/(gamma-1));

% output parameters to screen
fprintf('\n %20s %.15f %s \n','M2 =',M2,'');
fprintf('%20s %.15f %s \n','P2 =',P2,'Pa');
fprintf('%20s %.15f %s \n','T2 =',T2,'K');
fprintf('%20s %.15f %s \n','rho2 =',rho2,'kg/m^3');
fprintf('%20s %.15f %s \n','Pt2 =',Pt2,'Pa');

% clear all variables
clear;

% define plot parameters
gamma=1.4;
n=1000;
M1=linspace(1,3,n);

% calculate plot values
for i=1:n
    M2(i)=sqrt((M1(i)^2+2/(gamma-1))/(2*gamma/(gamma-1)*M1(i)^2-1));
    Pratio(i)=(1+gamma*M1(i)^2)/(1+gamma*M2(i)^2);
    Tratio(i)=(1+(gamma-1)/2*M1(i)^2)/(1+(gamma-1)/2*M2(i)^2);
    rhoratio(i)=M1(i)/M2(i)*sqrt((1+(gamma-1)/2*M2(i)^2)/(1+(gamma-1)/2*M1(i)^2));
    Ptratio(i)=(1+gamma*M1(i)^2)/(1+gamma*M2(i)^2)*((1+(gamma-1)/2*M2(i)^2)/(1+(gamma-1)/2*M1(i)^2))^(gamma/(gamma-1));
end

% construct plot
figure(1);
plot(M1,M2,'r');
hold on;
plot(M1,Pratio,'b');
plot(M1,rhoratio,'c');
plot(M1,Tratio,'m');
plot(M1,Ptratio,'g');
axis([1 3 0 5]);
legend('M2','P2/P1','rho2/rho1','T2/T1','Pt2/Pt1');
title('Normal Shock Jump Conditions');
xlabel('M1');