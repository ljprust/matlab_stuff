clear,clc;
h=linspace(0,15,1000);
q=.5*1.225*exp(h./(-6.7))*200^2;
plot(q,h);
title('Dynamic Pressure vs. Altitude');
xlabel('Dynamic Pressure (Pa)');
ylabel('Altitude (km)');