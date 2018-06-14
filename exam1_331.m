% Logan Prust - Aer E 331 - Exam 1 - Problem 2

clear;

wn=3.3;
z=0.01;
Iz=0.045;
wd=wn*sqrt(1-z^2);
psi0=5;
A=[0 1;-Iz*wn -2*z*Iz*wn];
B=[0;1];
C=[0 1];
D=0;
Gss=ss(A,B,C,D);
figure(1);
grid on;
hold on;
initial(Gss,[psi0,psi0]);

time=0:1:4000;
n=length(time);
for i=1:n
    psi(i)=psi0*exp(-wn^2*time(i))*(cos(wd*time(i))+wn^2/wd*sin(wd*time(i)));
end

plot(time,psi,'r');

clear;

Gp=tf([.73,1.43,.3,.24],[1,2.21,1.78,2.13,.07]);
figure(2);
grid on;
rlocus(Gp);

K=2.22;
W=feedback(K*Gp,1);

Gp0=3.4286;
W0=0.88388;

figure(3);
step(Gp/Gp0);
hold on;
step(W/W0);