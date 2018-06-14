% Logan Prust - Aer E 331 - Exam 3

clear;
format long;

%PROBLEM 1: [Etkin jet transport long. dynamics]

%PART A

%State vector: [u,w,q,th]  p.166 (6.2.1)

a1=[-0.006868, 0.013950,0.00000,-32.000];
a2=[-0.090550,-0.315100,773.980, 0.0000];
a3=[0.0001187,-0.001026,-0.4285, 0.0000];
a4=[0.0000000,0.0000000,1.00000, 0.0000];
A=[a1;a2;a3;a4];

%Control input: [de]  p.229 (7.6.5)

B=[-0.000187 ; -17.85 ; -1.158 ; 0];
C=eye(4); D=zeros(4,1);

[NUM,DEN]=ss2tf(A,B,C,D);

%PART B 

sys=ss(A,B,C,D);

figure(1);
impulse(sys*0.1,10);

%PART C

P=[-0.1 -0.1 -1+0.2i -1-0.2i];
K=acker(A,B,P);

%PART D

sys2=feedback(sys,K);

figure(2);
impulse(sys2*0.1,10);

%PROBLEM 2

%PART A

clear;

A=[-1 -9;1 0];
B=[2;0];
C=[0 1];
D=0;
P=[-5+5i -5-5i];

K=acker(A,B,P)

[NUM,DEN]=ss2tf(A,B,C,D);
G=tf(NUM,DEN);
s=tf('s');
Gc=K*[s;1];
What=feedback(G,Gc);
scale=50/2;
W=What*scale;

%PART B

figure(3);
rlocus(G*Gc);

%PART C

Wforhat=G*Gc/(1+G*Gc);
scale=450/369;
Wfor=scale*Wforhat;

figure(4);
hold on;
step(W);
step(Wfor);

%PART E

clear;

A=[-1 -9 0;1 0 0;0 1 0];
B=[2;0;0];
C=[0 1 0];
D=0;
P=[-5+5i -5-5i -5];

K=acker(A,B,P);

sys=ss(A,B,C,D);
s=tf('s');
Gc=K*[s;1;1/s];
W=feedback(sys,Gc);

figure(5);
step(W);

%PROBLEM 3

%PART C

z=tf('z',0.3142);
Gpz=(0.08287*z+0.07440)/(z^2+1.023*z+0.7304);
Gpzhat=0.4629*z/(z^2-1.023*z+0.7304);

figure(7);
hold on;
step(Gpz);
step(Gpzhat);

%PART D

clear;

Gps=tf(2,[1 1 9]);
Gpz1=c2d(Gps,0.3142);
Gpz2=c2d(Gps,0.6283);

figure(8);
hold on;
bode(Gps);
bode(Gpz1);
bode(Gpz2);

%PROBLEM 4

%PART C

T=0.0314;
s=tf('s');
z=tf('z',T);
Gcs=4.5*(s+4.5);
Gcz=4.5*(1/T+4.5-1/T/z);
Gcshat=d2c(Gcz,'tustin');

figure(9);
hold on;
bode(Gcs);
bode(Gcshat);

%PART D

Gps=tf(2,[1 1 9]);
Gs=Gps*Gcs;
Gshat=Gps*Gcshat;

figure(10);
hold on;
bode(Gs);
bode(Gshat);

%PART E

W=Gps/(1+Gs);
W2=feedback(Gps,Gcs);
What=Gps/(1+Gshat);
What2=feedback(Gps,Gcshat);

figure(11);
hold on;
step(W);
step(What);