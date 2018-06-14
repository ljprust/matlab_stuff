% Logan Prust - Aer E 331 - Homework 5

clear,clc;

% PROBLEM 1(c)

Ap=[-5 -3 -2;1 0 0;0 1 0];
Bp=[1;0;0];
Pp=[-4+3i -4-3i -4];

K = acker(Ap,Bp,Pp)

% PROBLEM 1(e)

Cp=[0 0 1];
Dp=0;

sysp=ss(Ap,Bp,Cp,Dp);
figure(1);
step(sysp);

figure(2);
hold on;
step(sysp);
step(2*sysp);

% PROBLEM 2(b)

A=[-0.76 -4.55;1 0];
eig(A)

roots([1,0.76,4.55])

% PROBLEM 2(d)

P=roots([1 1.54 1.06]);
B=[-4.6;0];
C=[0 1];
D=0;
acker(A,B,P)

sys=ss(A,B,C,D);

[NUM,DEN] = ss2tf(A,B,C,D);
G=tf(NUM,DEN)

scale=4.55/4.6;

figure(3);
hold on;
step(2*sysp);
step(scale*sys);

% PROBLEM 3(a)

clear;

A=[-.25 0 -1 .18;-16.02 -8.4 2.19 0;4.49 -.35 -.76 0;0 1 0 0];
B=[.07 0;23.16 -29.01;4.55 .22;0 0];
C=eye(4);
D=zeros(4,2);

eig(A)

[NUM,DEN]=ss2tf(A,B,C,D,2)

G=tf(NUM(4,:),DEN)

figure(4);
step(G);

% PROBLEM 4(c)

clear;

A=[-11.4 -14 0;1 0 0;0 1 0];
B=[500;0;0];

P=[-2,-2+1.5i,-2-1.5i];

K=place(A,B,P)

ABK=A-B*K;

% PROBLEM 4(d)

ks=1;
C=[0 0 ks];
D=0;

sys=ss(ABK,B,C,D);

figure(5);
step(sys);

ksnew=1/40;
Cnew=[0 0 ksnew];

sysnew=ss(ABK,B,Cnew,D);

figure(6);
step(sysnew);

% PROBLEM 5(b)

clear;

A=[-2 -25;1 0];
B=[50;0];
C=[0 1];
D=0;

sys=ss(A,B,C,D);
x0=[3 1];

figure(7);
initial(sys,x0);