clear,clc;
h=linspace(0,15,1000);
V=200*exp(h./6.7);
plot(V,h);
title('Airspeed vs. Altitude');
xlabel('Airspeed (m/s)');
ylabel('Altitude (km)');