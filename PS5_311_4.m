% Logan Prust - Aer E 311 - Homework 5 - Problem 4

clear,clc;

% set constants
alphadeg=linspace(-16,24,1000);
alpha=alphadeg*3.1415926536/180;
M=2.13;
k=0.007938;

% calculations
for i=1:1000
    Cl(i)=4*alpha(i)/sqrt(M^2-1);
    Cd(i)=4/sqrt(M^2-1)*(alpha(i)^2+k);
    LD(i)=alpha(i)/(alpha(i)^2+k);
end

% plots
figure(1);
plot(alphadeg,Cl);
xlabel('Angle of Attack (degrees)');
ylabel('Coefficient of Lift');
title('Coefficient of Lift vs Angle of Attack');
axis([-16 24 -0.8 1.2]);

figure(2);
plot(alphadeg,Cd);
xlabel('Angle of Attack (degrees)');
ylabel('Coefficient of Drag');
title('Coefficient of Drag vs Angle of Attack');
axis([-16 24 0 0.4]);

figure(3);
plot(alphadeg,LD);
xlabel('Angle of Attack (degrees)');
ylabel('Lift Over Drag');
title('Lift Over Drag vs Angle of Attack');
axis([-16 24 -8 8]);