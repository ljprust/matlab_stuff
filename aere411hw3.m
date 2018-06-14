% Logan Prust - Aer E 411 - Homework 3 - Problem 2

clear;

% define givens
ma = 800; % kg/s
Tt1 = 500; % K
hf = 43000*1000; % J/kg
mf = 8; % kg/s
R = 287; % J/kg/K
gamma = 1.4;

% calculate some stuff
Cp = gamma*R/(gamma-1);
k = mf*hf/ma/Cp/Tt1 + 1;

% initialize vectors
n=1000;
M1 = linspace(0,0.4,n);
M2 = zeros(1,n);
pib = zeros(1,n);

% solve
for i=1:n
    A = (gamma-1)/2 + 0.5*gamma*(gamma-1)*M1(i)^2;
    B = 1 + gamma*M1(i)^2 - k*gamma*M1(i)^2;
    C = -k*M1(i)^2;
    M2(i) = sqrt((-B + sqrt(B^2 - 4*A*C))/2/A);
    pib(i) = (1+gamma*M1(i)^2)/(1+gamma*M2(i)^2)*((1+(gamma-1)/2*M2(i)^2)/(1+(gamma-1)/2*M1(i)^2))^(gamma/(gamma-1));
end

% plot output
figure(1);
hold on;
plot(M1,M2,'r');
plot(M1,pib);
title('Pressure Ratio and Downstream Mach Number vs Upstream Mach Number');
xlabel('Upstream Mach Number M1');
legend('M2','pi');
axis([0 0.4 0 1.3]);