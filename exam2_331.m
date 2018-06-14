% Logan Prust - Aer E 331 - Exam 2

clear,clc;

% PROBLEM 1(a)

Gzeros=roots([.7 1.4 .3 .2])
Gpoles=roots([1 2.2 1.8 2.1 .4])

% PROBLEM 1(b)

s=tf('s');
Gp=tf([.7 1.4 .3 .2],[1 2.2 1.8 2.1 .4]);
Gpnew=0.727*(s^2+2*s*0.0741+0.0741^2+0.386^2)/(s^2+2*s*0.0984+0.0984^2+1)/(s+0.2225);

figure(1);
hold on;
impulse(Gp);
impulse(Gpnew);

% PROBLEM 1(c)

figure(2);
hold on;
bode(Gp);
bode(Gpnew);

% PROBLEM 1(d)

figure(3);
step(Gp*10);

% PROBLEM 2(c)

s=tf('s');
Gm=120/s/(s+1);
Gc=(s+1)/(s+5);
G=Gm*Gc;
figure(4);
rlocus(G);

gamma=2.375;
t=0:0.01:5;
u=gamma*t;
figure(5);
lsim(G,u,t);

% PROBLEM 2(d)

data=load('ICBM.txt');
t=data(1,:);
y=data(2,:);
figure(6);
lsim(G,y,t);

% PROBLEM 3(a)

s=tf('s');
K=1/30;
Gm=120/s/(s+1);
figure(7);
bode(K*Gm);

% PROBLEM 3(b)

Gc=3.031/(3.031-s);
figure(8);
bode(Gc*K*Gm);

% PROBLEM 3(c)

figure(9);
lsim(Gc*K*Gm,y,t);

% PROBLEM 3(d)

figure(10);
hold on;
step(G,1);
step(Gc*K*Gm,1);
legend('Problem 2 Controller','Problem 3 Controller');

% PROBLEM 4(b)

delta=0.1;
s=tf('s');
Gc=tf([.1174,1],[.0155,1]);
Gp=100*exp(-s*2*delta)/s/(.5*s+1);

figure(11);
bode(Gc*Gp);