% Logan Prust - Aer E 411 - Homework 5 - Problem 3

clear,clc;

% define constants
pic = 10;
gamma = 1.4;
T0 = 210;
Cp = 1005;
h = 4.42E7;
taulambda = 7;
taulambdaAB = 8;

% calculate stuff
R = Cp*(gamma-1)/gamma;
a0 = sqrt(gamma*R*T0);
tauc = pic^((gamma-1)/gamma);

n=1000;
M0 = linspace(0,3.5,n);

% carry out other calculations
for i=1:n
    taur(i) = 1 + (gamma-1)/2*M0(i)^2;
    taut(i) = 1 - taur(i)/taulambda*(tauc-1);
    spcthrust(i) = a0*(sqrt(2/(gamma-1)*taulambdaAB/taut(i)/tauc/taur(i)*(taut(i)*tauc*taur(i)-1)) - M0(i));
    ftotal(i) = Cp*T0/h*(taulambdaAB - taur(i));
    S(i) = ftotal(i)/spcthrust(i);
end

% plot results
figure(1);
plot(M0,spcthrust);
xlabel('Mach Number');
ylabel('Specific Thrust (m/s)');
title('Specific Thrust vs Mach Number');

figure(2);
plot(M0,S);
xlabel('Mach Number');
ylabel('Specific Fuel Consumption (s/m)');
title('Specific Fuel Consumption vs Mach Number');