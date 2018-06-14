% Logan Prust - Phys 322L - NMR Lab

clear;

format long;

k=89.66; %uT/A

x=[15 28 45 57 76 89 105 122 135 0 -16 -32 -44 -62 -76 -93 -105]; % mA
y=[2025.5 2071.5 2137.5 2178 2254 2300 2386.5 2427 2472 1975 1913.5 1851 1807.5 1740.5 1689.5 1623 1578]; % Hz

n=length(x);

for i=1:n
    X(i,1)=x(i);
    X(i,2)=1;
end

beta=inv(transpose(X)*X)*transpose(X)*y';

m=beta(1)

xbar=mean(x);
ybar=mean(y);

ssxx=0;
ssyy=0;
ssxy=0;
for i=1:n
    ssxx=ssxx+(x(i)-xbar)^2;
    ssyy=ssyy+(y(i)-ybar)^2;
    ssxy=ssxy+(x(i)-xbar)*(y(i)-ybar);
end

s=sqrt((ssyy-beta(1)*ssxy)/(n-2));

SEa=s*sqrt(1/n+xbar^2/ssxx);

SEb=s/sqrt(ssxx);

sigmam=SEb

c=m/k*1000

sigmac=sigmam/k*1000

f0=beta(2)

sigmaf0=SEa

error=0.25*ones(1,n);
y1=beta(1)*(-150)+beta(2);
y2=beta(1)*(150)+beta(2);
figure(2);
hold on;
errorbar(x,y,error,'*');
plot([-150 150],[y1 y2],'r');
title('Frequency vs Current');
xlabel('Current (mA)');
ylabel('Frequency (Hz)');

polartime=0.5:0.5:8;
polarwater=10^(-5)*[.4 1.4 2.7 3.7 6.5 7.4 8 9 11 13 11 13 15 13 16 14];
polarPFPE=10^(-6)*[20 22 14 4 4.4 4.7 4.4 4 4.2 4.5 7.6 7.8 8.2 8 7.8 7];
polarwatererror=10^(-5)*[.005 .05 .05 .05 .05 .05 .05 .05 .5 .5 .5 .5 .5 .5 .5 .5];
polarPFPEerror=10^(-6)*[.5 .5 .5 .05 .05 .05 .05 .05 .05 .05 .05 .05 .05 .05 .05 .05];

figure(1);
hold on;
errorbar(polartime,polarwater,polarwatererror,'*');
vline=get(gca,'ylim');
plot([5 5],vline,'r');
title('Polarization Convergence for Tap Water');
xlabel('Polarization Time (sec)');
ylabel('Amplitude');

figure(3);
hold on;
errorbar(polartime,polarPFPE,polarPFPEerror,'*');
vline2=get(gca,'ylim');
plot([5.5 5.5],vline2,'r');
title('Polarization Convergence for PFPE');
xlabel('Polarization Time (sec)');
ylabel('Amplitude');