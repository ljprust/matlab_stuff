% Logan Prust Standard Atmosphere Model Program (Aer E 161 Project 1)
clear,clc;
h=0:1:15; % set up geometric altitude matrix
I=dlmread('AtmosphericInfo.txt',',',[4,0,4,0]); % get value at which isothermal layer starts
T0=dlmread('AtmosphericInfo.txt',',',[7,0,7,0]); % get temperature at sea level
p0=dlmread('AtmosphericInfo.txt',',',[10,0,10,0]); % get pressure at sea level
roh0=dlmread('AtmosphericInfo.txt',',',[13,0,13,0]); % get density at sea level
g0=dlmread('AtmosphericInfo.txt',',',[16,0,16,0]); % get acceleration of gravity at sea level
R=dlmread('AtmosphericInfo.txt',',',[19,0,19,0]); % get gas constant
a=dlmread('AtmosphericInfo.txt',',',[22,0,22,0]); % get lapse rate
rE=dlmread('AtmosphericInfo.txt',',',[25,0,25,0]); % get radius of Earth
hG=h(1:(I+1)); % get geometric altitude values for gradient layer
hI=h((I+1):(15+1)); % get geometric altitude values for isothermal layer
TG=SAM_G_T(T0,a,hG); % calculate temperature values for gradient layer
pG=SAM_G_p(p0,TG,T0,g0,a,R); % calculate pressure values for gradient layer
rohG=SAM_G_roh(roh0,TG,T0,g0,a,R); % calculate density values for gradient layer
TI=SAM_I_T(TG); % calculate temperature values for isothermal layer
pI=SAM_I_p(pG,g0,R,TG,hI,I); % calculate pressure values for isothermal layer
rohI=SAM_I_roh(rohG,g0,R,TG,hI,I); % calculate density values for isothermal layer
T=[TG,TI]; % assemble temperature matrix
p=[pG,pI(2:length(pI))]; % assemble pressure matrix
roh=[rohG,rohI(2:length(rohI))]; % assemble density matrix
theta=T/T0; % calculate theta values
delta=p/p0; % calculate delta values
sigma=roh/roh0; % calculate sigma values
geopotalt=rE*h./(rE+h); % calculate geopotential altitude values
SAT=[h',geopotalt',T',p'/1000,roh',theta',delta',sigma']; % create Standard Atmosphere Table
fprintf('%10s %8s %8s %9s %13s %5s %9s %9s\n','  h(g) (km)','h (km)','T (K)','p (kPa)','roh (kg/m^3)','theta','delta','sigma'); % output column headings to screen
disp(' '); % leave a line empty to make it prettier
disp(SAT); % output Standard Atmosphere Table to screen
% make plots of all variables with respect to geopotential altitude
figure(1);
plot(T,geopotalt);
title('Temperature vs. Geopotential Altitude');
xlabel('Temperature (K)');
ylabel('Geopotential Altitude (km)');
figure(2);
plot(roh,geopotalt);
title('Density vs. Geopotential Altitude');
xlabel('Density (kg/m^3)');
ylabel('Geopotential Altitude');
figure(3);
plot(p,geopotalt);
title('Pressure vs. Geopotential Altitude');
xlabel('Pressure (Pa)');
ylabel('Geopotential Altitude');
figure(4);
plot(theta,geopotalt);
title('theta vs. Geopotential Altitude');
xlabel('theta (dimensionless)');
ylabel('Geopotential Altitude');
figure(5);
plot(sigma,geopotalt);
title('sigma vs. Geopotential Altitude');
xlabel('sigma (dimensionless)');
ylabel('Geopotential Altitude');
figure(6);
plot(delta,geopotalt);
title('delta vs. Geopotential Altitude');
xlabel('delta (dimensionless)');
ylabel('Geopotential Altitude');