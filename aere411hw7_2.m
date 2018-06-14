% Logan Prust - Aer E 411 - Homework 7 - Problem 2

clear,clc;

% define constants
M3 = 0.1;
gammac = 1.4;
gammat = 1.28;
CD = 1.5;

% initialize total temperature ratio
n = 1000;
Tt = linspace(1.5,3.5,n);

% carry out calculations
for i=1:n
    chi = gammac/gammat*M3^2*(1+(gammac-1)/2*M3^2)/(1+gammac*M3^2*(1-CD/2))^2*Tt(i);
    M4(i) = sqrt(2*chi/(1-2*gammat*chi+sqrt(1-2*(gammat+1)*chi)));
    P = (1+gammac*M3^2*(1-CD/2))/(1+gammat*M4(i)^2);
    pib(i) = P*(1+(gammat-1)/2*M4(i)^2)^(gammat/(gammat-1))/(1+(gammac-1)/2*M3^2)^(gammac/(gammac-1));
end

% plot results
figure(1);
plot(Tt,M4);
xlabel('Total Temperature Ratio');
ylabel('Exit Mach Number M4');
title('Exit Mach Number vs Total Temperature Ratio');

figure(2);
plot(Tt,pib);
xlabel('Total Temperature Ratio');
ylabel('Total Pressure Ratio');
title('Total Pressure Ratio vs Total Temperature Ratio');