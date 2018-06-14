% Logan Prust - Aer E 411 - Homework 4 - Problem 3

clear,clc;

% define givens
gamma = 1.4;
T0 = 222.2;
h = 4.4194E7;
Cp = 1004.5;
taulambda = 6.8;

% calculate some stuff
R = (gamma-1)/gamma*Cp;
a0 = sqrt(gamma*R*T0);
M0max = sqrt(2/(gamma-1)*(sqrt(taulambda)-1));
M0op = sqrt(2/(gamma-1)*(taulambda-1));

% initalize M0 vector
n = 1000;
M0 = linspace(0,M0op,n);

% calculate other stuff
for i = 1:n
    spcthrust(i) = a0*(sqrt(2/(gamma-1)*(sqrt(taulambda)-1)^2+M0(i)^2)-M0(i));
    taur(i) = 1+(gamma-1)/2*M0(i)^2;
    pic(i) = (sqrt(taulambda)/taur(i))^(gamma/(gamma-1));
    f(i) = Cp*T0/h*(taulambda - taur(i));
    S(i) = f(i)/spcthrust(i);
end

% plot results
figure(1);
plot(M0,spcthrust);
xlabel('Mach Number');
ylabel('Specific Thrust (N*s/kg)');
title('Specific Thrust versus Mach Number');

figure(2);
hold on;
plot(M0,pic);
plot([0 6],[1 1]);
plot([M0max M0max],[0 30]);
xlabel('Mach Number');
ylabel('Compressor Ratio');
title('Compressor Ratio versus Mach Number');

figure(3);
plot(M0,S);
xlabel('Mach Number');
ylabel('Specific Fuel Consumption');
title('Specific Fuel Consumption versus Mach Number');