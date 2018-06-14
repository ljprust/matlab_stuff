% Logan Prust - Aer E 411 - Homework 4 - Problem 1

clear,clc;

% define givens
M0 = 2.7;
gamma = 1.4;
T0 = 222.2;
h = 4.4194E7;
Cp = 1004.5;
taulambda = 6.8;

% calculate some things
taur = 1+(gamma-1)/2*M0^2;
R = (gamma-1)/gamma*Cp;
a0 = sqrt(gamma*R*T0);

% initialize pi_c vector
n=1000;
pic = linspace(0,40,n);

% calculate other things
for i=1:n
    tauc(i) = pic(i)^((gamma-1)/gamma);
    f(i) = Cp*T0/h*(taulambda - tauc(i)*taur);
    taut(i) = 1-taur/taulambda*(tauc(i) - 1);
    spcthrust(i) = a0*(sqrt(2/(gamma-1)*taulambda/tauc(i)/taur*(taut(i)*tauc(i)*taur-1))-M0);
    S(i) = f(i)/spcthrust(i);
end

% plot the results
figure(1);
plot(pic,spcthrust);
xlabel('Pressure Ratio');
ylabel('Specific Thrust (N*s/kg)');
title('Specific Thrust versus Pressure Ratio');
axis([0 40 0 600]);

figure(2);
plot(pic,S);
xlabel('Pressure Ratio');
ylabel('Specific Fuel Consumption');
title('Specific Fuel Consumption versus Pressure Ratio');
axis([0 40 0 0.0001]);