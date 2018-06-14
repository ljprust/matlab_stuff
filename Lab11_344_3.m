clear all
clc
x = [-1.50,-.30,-.18,0.00,.15,.30,.45,.60,.75,.90,1.05,1.20,1.35,1.45];
Me = [0.6779 0.9025 0.9304 1.000 1.243 1.348 1.431 0.6699 0.6254 0.5936 0.5675 0.5521 0.5455 0.5429];

Pe = 97924.4; %atmosphere pressure = exit pressure
y = 1.4; %gamma 

M2 =( 1 + (((y-1)/2).*Me.^2))./((y.*Me.^2)-((y-1)/2));


Pt = ((1 + ((y-1)/2).*(M2.^2)).^(y/(y-1))).*Pe; %total pressure behind shock
P = Pt.*((1 + (((y-1)/2).*(Me.^2))).^(-y/(y-1))); %Pressure behind shock

P = P./1000;

figure(1)
plot(x,P)
axis([-4 2 0 6e4]);
title('1st Critical')
xlabel('Nozzle Location (in)')
ylabel(' P (Pa)') 
figure(2)
plot(x,Pt)
axis([-4 2 0 6e7]);
title('1st Critical')
xlabel('Nozzle Location (in)')
ylabel(' Pt (Pa)') 
