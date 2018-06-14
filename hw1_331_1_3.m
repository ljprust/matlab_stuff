% Logan Prust - Aer E 331 - Homework 1 - Problem 3

clear,clc;

s=tf('s');
B=(.0997*s-.4822)/(s^2+.1455*s+.4822);
T1=feedback(ss(B),0);
figure(1);
step(T1,100);

R=(.5101*s)/(s^2+.1455*s+.4822);
T2=feedback(ss(R),0);
figure(2);
step(T2,100);

PSI=(.5101)/(s^2+.1455*s+.4822);
T3=feedback(ss(PSI),0);
figure(3);
step(T3,100);

figure(4);
hold on;
step(T3,100);
step(-T1,100);

H=-PSI/B;
T4=feedback(ss(H),0);
figure(5);
bode(T4);