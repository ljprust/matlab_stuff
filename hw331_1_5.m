% Logan Prust - Aer E 331 - Homework 1 - Problem 3

clear,clc;

s=tf('s');
Gp=10/(s^2+0.2*s+100);
Kp=190;
Kd=8.924;
Gc=Kp+Kd*s;
W=Gc*Gp/(1+Gc*Gp);

TGp=feedback(ss(Gp),0);
figure(1);
step(TGp,100);

TW=feedback(ss(W),0);
figure(2);
step(TW,0.5);