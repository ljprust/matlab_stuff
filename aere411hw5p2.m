% Logan Prust - Aer E 411 - Homework 5 - Problem 2

clear,clc;

% define constants
gamma = 1.4;
T0 = 220;
Cp = 1005;
M0 = 2.5;
h = 4.42E7;

% calculate stuff
R = Cp*(gamma-1)/gamma;
a0 = sqrt(gamma*R*T0);
taur = 1+(gamma-1)/2*M0^2;

% Part (b)

taulambda = 7;
n=1000;
taulambdaAB = linspace(7,9,n);

% carry out all calculations
for i=1:n
    spcthrust(i) = a0*(sqrt(2*taulambdaAB(i)/(gamma-1)*(1-4*taulambda/(taulambda+taur)^2)) - M0);
    ftotal(i) = Cp*T0/h*(taulambdaAB(i) - taur);
    S(i) = ftotal(i)/spcthrust(i);
end

% plot results
figure(1);
plot(taulambdaAB,spcthrust);
xlabel('tau lambda AB');
ylabel('Specific Thrust (m/s)');
title('Specific Thrust vs Temperature Ratio of the Afterburner');

figure(2);
plot(taulambdaAB,S);
xlabel('tau lambda AB');
ylabel('Specific Fuel Consumption (s/m)');
title('Specific Fuel Consumption vs Temperature Ratio of the Afterburner');

% Part (c)

% clear necessary variables
clear taulambda taulambdaAB spcthrust S ftotal;

taulambdaAB = 8;
taulambda = linspace(6.5,7.5,n);

% calculate ftotal
ftotal = Cp*T0/h*(taulambdaAB - taur);

% carry out other calculations
for i=1:n
    spcthrust(i) = a0*(sqrt(2*taulambdaAB/(gamma-1)*(1-4*taulambda(i)/(taulambda(i)+taur)^2)) - M0);
    S(i) = ftotal/spcthrust(i);
end

% plot results
figure(3);
plot(taulambda,spcthrust);
xlabel('tau lambda');
ylabel('Specific Thrust (m/s)');
title('Specific Thrust vs Temperature Ratio of tau lambda');

figure(4);
plot(taulambda,S);
xlabel('tau lambda');
ylabel('Specific Fuel Consumption (s/m)');
title('Specific Fuel Consumption vs Temperature Ratio of tau lambda');