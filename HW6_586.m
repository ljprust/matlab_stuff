% Logan Prust - Astro 586 - Homework 6 - Problem 2

clear,clc;

R0=1;
vc=300;
alpha1=-0.2;
alpha2=0;
alpha3=0.2;
n=100;
R=linspace(0,30,n);
c=5;
G=(6.673E-11)/1000^2/(3.0857E19); % km^2/s^2*kpc/kg
Rcr=10; % kpc

pattern1=vc/(R0+Rcr^(1+0.2)); % km/s/kpc
pattern2=vc/(R0+Rcr^(1+0));
pattern3=vc/(R0+Rcr^(1-0.2));

for i=1:n
    v1(i)=vc*R(i)/(R0+R(i)^(1-alpha1));
    v2(i)=vc*R(i)/(R0+R(i)^(1-alpha2));
    v3(i)=vc*R(i)/(R0+R(i)^(1-alpha3));
    alpha=alpha1;
    Emin=c/pi/G*vc*sqrt(4/(R0+R(i)^(1-alpha))^2-2*(1-alpha)*R(i)^(-alpha)/(R0+R(i)^(1-alpha))^3);
    m=2;
    a1m2p(i)=vc/(R0+R(i)^(1-alpha))+vc/m*sqrt(4/(R0+R(i)^(1-alpha))^2-2*(1-alpha)*R(i)^(-alpha)/(R0+R(i)^(1-alpha))^3);
    a1m2m(i)=vc/(R0+R(i)^(1-alpha))-vc/m*sqrt(4/(R0+R(i)^(1-alpha))^2-2*(1-alpha)*R(i)^(-alpha)/(R0+R(i)^(1-alpha))^3);
    m=3;
    a1m3p(i)=vc/(R0+R(i)^(1-alpha))+vc/m*sqrt(4/(R0+R(i)^(1-alpha))^2-2*(1-alpha)*R(i)^(-alpha)/(R0+R(i)^(1-alpha))^3);
    a1m3m(i)=vc/(R0+R(i)^(1-alpha))-vc/m*sqrt(4/(R0+R(i)^(1-alpha))^2-2*(1-alpha)*R(i)^(-alpha)/(R0+R(i)^(1-alpha))^3);
    alpha=alpha2;
    m=2;
    a2m2p(i)=vc/(R0+R(i)^(1-alpha))+vc/m*sqrt(4/(R0+R(i)^(1-alpha))^2-2*(1-alpha)*R(i)^(-alpha)/(R0+R(i)^(1-alpha))^3);
    a2m2m(i)=vc/(R0+R(i)^(1-alpha))-vc/m*sqrt(4/(R0+R(i)^(1-alpha))^2-2*(1-alpha)*R(i)^(-alpha)/(R0+R(i)^(1-alpha))^3);
    m=3;
    a2m3p(i)=vc/(R0+R(i)^(1-alpha))+vc/m*sqrt(4/(R0+R(i)^(1-alpha))^2-2*(1-alpha)*R(i)^(-alpha)/(R0+R(i)^(1-alpha))^3);
    a2m3m(i)=vc/(R0+R(i)^(1-alpha))-vc/m*sqrt(4/(R0+R(i)^(1-alpha))^2-2*(1-alpha)*R(i)^(-alpha)/(R0+R(i)^(1-alpha))^3);
    alpha=alpha3;
    m=2;
    a3m2p(i)=vc/(R0+R(i)^(1-alpha))+vc/m*sqrt(4/(R0+R(i)^(1-alpha))^2-2*(1-alpha)*R(i)^(-alpha)/(R0+R(i)^(1-alpha))^3);
    a3m2m(i)=vc/(R0+R(i)^(1-alpha))-vc/m*sqrt(4/(R0+R(i)^(1-alpha))^2-2*(1-alpha)*R(i)^(-alpha)/(R0+R(i)^(1-alpha))^3);
    m=3;
    a3m3p(i)=vc/(R0+R(i)^(1-alpha))+vc/m*sqrt(4/(R0+R(i)^(1-alpha))^2-2*(1-alpha)*R(i)^(-alpha)/(R0+R(i)^(1-alpha))^3);
    a3m3m(i)=vc/(R0+R(i)^(1-alpha))-vc/m*sqrt(4/(R0+R(i)^(1-alpha))^2-2*(1-alpha)*R(i)^(-alpha)/(R0+R(i)^(1-alpha))^3);
end

for i=1:n
    alpha=alpha1;
    kcrit1(i)=1/c*vc*sqrt(4/(R0+R(i)^(1-alpha))^2-2*(1-alpha)*R(i)^(-alpha)/(R0+R(i)^(1-alpha))^3);
    kcritbar1(i)=kcrit1(i)*R(i)/2;
    x1(i)=kcritbar1(i)/2;
    alpha=alpha2;
    kcrit2(i)=1/c*vc*sqrt(4/(R0+R(i)^(1-alpha))^2-2*(1-alpha)*R(i)^(-alpha)/(R0+R(i)^(1-alpha))^3);
    kcritbar2(i)=kcrit2(i)*R(i)/2;
    x2(i)=kcritbar2(i)/2;
    alpha=alpha3;
    kcrit3(i)=1/c*vc*sqrt(4/(R0+R(i)^(1-alpha))^2-2*(1-alpha)*R(i)^(-alpha)/(R0+R(i)^(1-alpha))^3);
    kcritbar3(i)=kcrit3(i)*R(i)/2;
    x3(i)=kcritbar3(i)/2;
end

figure(1);
hold on;
plot(R,v1,'r');
plot(R,v2,'b');
plot(R,v3,'g');
title('Rotation Curves');
xlabel('Radius (kpc)');
ylabel('Velocity (km/s)');
legend('alpha = -0.2','alpha = 0.0','alpha = 0.2');

figure(2);
hold on;
plot(R,a1m2p,'r');
plot(R,a1m2m,'b');
plot(R,a1m3p,'g');
plot(R,a1m3m,'k');
plot(R,pattern1,'-');
title('Alpha = -0.2');
xlabel('Radius (kpc)');
ylabel('Frequency (km/s/kpc)');
legend('m=2, +','m=2, -','m=3, +','m=3, -');

figure(3);
hold on;
plot(R,a2m2p,'r');
plot(R,a2m2m,'b');
plot(R,a2m3p,'g');
plot(R,a2m3m,'k');
plot(R,pattern2,'-');
title('Alpha = 0.0');
xlabel('Radius (kpc)');
ylabel('Frequency (km/s/kpc)');
legend('m=2, +','m=2, -','m=3, +','m=3, -');

figure(4);
hold on;
plot(R,a3m2p,'r');
plot(R,a3m2m,'b');
plot(R,a3m3p,'g');
plot(R,a3m3m,'k');
plot(R,pattern3,'-');
title('Alpha = 0.2');
xlabel('Radius (kpc)');
ylabel('Frequency (km/s/kpc)');
legend('m=2, +','m=2, -','m=3, +','m=3, -');

figure(5);
hold on;
plot(R,kcrit1,'r');
plot(R,kcrit2,'b');
plot(R,kcrit3,'g');
title('Critical Wavenumber');
xlabel('Radius (kpc)');
ylabel('Critical Wavenumber (1/kpc)');
legend('alpha = -0.2','alpha = 0.0','alpha = 0.2');

figure(6);
hold on;
plot(R,x1,'r');
plot(R,x2,'b');
plot(R,x3,'g');
title('Bar Instability');
xlabel('Radius (kpc)');
ylabel('Instability Parameter');
legend('alpha = -0.2','alpha = 0.0','alpha = 0.2');