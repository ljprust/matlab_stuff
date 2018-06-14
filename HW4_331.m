% Logan Prust - Aer E 331 - Homework 4

% 1a

clear,clc;

frequency=[20 50 100 200 500 1000 1500 2000 3000 4000 5000 10000 15000 20000];
logfreq=log10(frequency);
phon100=[128 113 106 101 100 100 102 100 97 98 104 114 109 130];
phon60=[109 89 78 70 62 60 63 60 57 58 67 73 68 77];
phon2=[73 42 24 13 4 2 3 -2 -7 -6 3 13 12 20];

phon100=100-phon100;
phon60=60-phon60;
phon2=-phon2;

figure(1);
hold on;
plot(logfreq,phon100,'b');
plot(logfreq,phon60,'g');
plot(logfreq,phon2,'r');
grid on;
legend('100 phon','60 phon','threshold');

% 1a

clear,clc;

f=[20 50 100 200 500 1000 1500 2000 3000 4000 5000 10000 15000 20000]; %(Hz)
c1=100; c2=60; c3=2;
% The values in each of the three contours:
L1=[128 113 106 101 100 100 102 100 97 98 104 114 109 130];
L2=[109 89 78 70 62 60 63 60 57 58 67 73 68 77];
L3=[73 42 24 13 4 2 3 -2 -7 -6 3 13 12 20];
A1=c1-L1; A2=c2-L2; A3=c3-L3;
figure(1)
semilogx(f,[A1;A2;A3])
grid
legend('100 phon','60 phon','threshold','Location','SouthEast')
ylabel('dB')
xlabel('Frequency (Hz)')
title('Ear FRF Magnitudes at Three Phon Levels')

% 2

f1=[40,60,100,200,250,300,400,500,600,800,1000,1200,2000,2400,3000];
dB1=[50,35,24,12, 10,  9,  7,  6, 5.5, 4,  4,   4,   3.5, -1,  -8]; 
d1=[12 16 24 24 25 25 25 25 25 25 25 25 25 25 25];
%NOTE: Because the first 4 frequencies include only 12 (not 25) counts,
% they need to be multiplied by a factor of ~2.
f2=[3.5,4,4.5, 5, 5.5, 6, 7, 8, 9, 10, 11, 12, 14, 18]*10^3;
dB2=[-10,-8,-3, 0,  4,  7, 10,17,19, 17,12, 11, 15,22];
d2=[25 25 24 22 19 17 17 16 17 17 19 22 16 14];
d=[d1,d2]; feq=[f1,f2]; dB=[dB1,dB2];
d35=25;
eqL=3*(d35 - d)-4;
figure(7)
semilogx(feq,[eqL;dB])
grid on;
legend('Logan Threshold','Typical Threshold','Location','NorthEast')
ylabel('dB')
xlabel('Frequency (Hz)')
title('Logan & Typical Hearing Contours')

% 3a

clear,clc;

s=tf('s');
Gc=0.2*(s+0.05);
Gp=1/s^2;
G=Gp*Gc;
figure(2);
bode(G);

% 3b

Wc=feedback(Gc*Gp,1);
Wd=feedback(Gp,Gc);
figure(3);
bode(Wd);

% 4a

clear,clc;

w=[1 2 5 10 20 50 100];
Gmag=1./sqrt(w.^2+100)
Gphase=-atand(w/10)

% 5a

clear,clc;

s=tf('s');
Gp=10/s/(s+2);
Gc=4.13*(s+2.75)/(s+9.55);
figure(4);
bode(Gp);
figure(5);
bode(Gc);
figure(6);
bode(Gp*Gc);