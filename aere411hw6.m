% Logan Prust - Aer E 411 - Homework 6 - Problem 3b

clear,clc;

% define constants
Tt4a = 1500;
A8 = 0.1;
h = 4.4E7;
gamma = 1.4;
R = 287;
Cp = 1004.9;
M0 = 0.85;
T0 = 220;
P0 = 25000;
pic = 40;
picp = 2;

% calculate stuff
taur = 1+(gamma-1)/2*M0^2;
tauc = pic^((gamma-1)/gamma);
taucp = picp^((gamma-1)/gamma);
taulambda = Tt4a/T0;
tautc = 1-(tauc-1)*taur/taulambda;
f = Cp*T0/h*(taulambda-tauc*taur);
a0 = sqrt(gamma*R*T0);

% initiate alpha values
n = 1000;
alpha = linspace(0,20,n);

% carry out calculations
for i=1:n
    tautf = 1 - alpha(i)*(taucp-1)*taur/tautc/taulambda;
    spcthrust = a0*M0/(1+alpha(i))*((1/M0*sqrt(2/(gamma-1)*taulambda/taur/tauc*(tautc*tautf*tauc*taur-1))-1)+alpha(i)*(1/M0*sqrt(2/(gamma-1)*(taucp*taur-1))-1));
    S(i) = f/(1+alpha(i))/spcthrust;
end

% plot results
plot(alpha,S);
title('Specific Fuel Consumption vs Bypass Ratio');
xlabel('Bypass Ratio');
ylabel('Specific Fuel Consumption (s/m)');