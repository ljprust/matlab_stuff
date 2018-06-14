% Logan Prust - Aer E 311 - Homework 2 - Problem 6 - Part (a)

clear,clc;

gamma=1.4;

% get input parameters
M1=input('Upstream Mach Number: ');
sigma=input('Shock Angle (degrees): ');

% calculate output parameters
delta=atand(2*cotd(sigma)*(M1^2*(sind(sigma))^2-1)/(M1^2*(gamma+cosd(2*sigma))+2));
M1n=M1*sind(sigma);
M2n=(M1n^2+2/(gamma-1))/(2*gamma/(gamma-1)*M1n^2-1);
M2=M2n/sind(sigma-delta);

% print output parameters
fprintf('\n %20s %.15f %s \n','M2 =',M2,'');
fprintf('%20s %.15f %s \n','delta =',delta,'degrees');

clear;

% set parameters for plot
gamma=1.4;
n=1000;
M1=[1.5,2];
mu=[asind(1/1.5),asind(1/2)];

% create plot
for i=1:2
    sigma=linspace(mu(i),90,n);
    for j=1:n
        delta(j)=atand(2*cotd(sigma(j))*(M1(i)^2*(sind(sigma(j)))^2-1)/(M1(i)^2*(gamma+cosd(2*sigma(j)))+2));
    end
    if i==1
        plot(delta,sigma,'b');
    else
        plot(delta,sigma,'g');
    end
    hold on;
end
axis([0 25 20 100]);
legend('M = 1.5','M = 2');
title('Shock Angle vs Turning Angle');
xlabel('Turning Angle (degrees)');
ylabel('Shock Angle (degrees)');