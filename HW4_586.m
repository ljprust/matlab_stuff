% Logan Prust - Astro 586 - Homework 4

clear,clc;

psi=10;
psi2c=10;
R=0.17*10/psi;
tau1a=10/0.9;
tau1c=1/0.9;
tau2b=3;
n=100;
t=linspace(0,10,n);
Zsun=0.02;
y=Zsun;
Z0=0.1*Zsun;
Zf=Zsun;
M=100;

for i=1:n
    mu1a(i)=1-t(i)/tau1a;
    Z1a(i)=y*log(1/mu1a(i));
    Zs1a(i)=y*(1-mu1a(i)/(1-mu1a(i))*log(1/mu1a(i)));
    Z1b(i)=y*log(1/mu1a(i))+Z0;
    Zs1b(i)=y-mu1a(i)/(1-mu1a(i))*Z1b(i);
    mu1c(i)=1/(1+t(i)/tau1c);
    nu1c(i)=t(i)/tau1c;
    Z1c(i)=y*(1-exp(-nu1c(i)));
    Zs1c(i)=y-mu1c(i)/(1-mu1c(i))*Z1c(i);
    Z1d(i)=y+Zf-(y+Zf)*exp(-nu1c(i));
    Zs1d(i)=y-mu1c(i)/(1-mu1c(i))*Z1d(i);
    Z2a(i)=Z1a(i);
    Zs2a(i)=Zs1a(i);
    mu2b(i)=0.3111*tau2b*(exp(-t(i)/tau2b)-1)+1;
    if t(i) <= 1
        mu2c(i)=1;
    elseif t(i) <= 2
        mu2c(i)=1-0.3*(t(i)-1);
    elseif t(i) <= 3
        mu2c(i)=1-0.3*1;
    elseif t(i) <= 4
        mu2c(i)=1-0.3-0.3*(t(i)-3);
    elseif t(i) <= 5
        mu2c(i)=1-0.3-0.3;
    elseif t(i) <= 6
        mu2c(i)=1-0.3-0.3-0.3*(t(i)-5);
    else
        mu2c(i)=0.1;
    end
    if t(i) < 3
        mu2d(i)=-2/15*(t(i)^3-0.25*t(i)^4)+1;
    else
        mu2d(i)=0.1;
    end
    Z2b(i)=y*log(1/mu2b(i));
    Zs2b(i)=y-mu2b(i)/(1-mu2b(i))*Z2b(i);
    Z2c(i)=y*log(1/mu2c(i));
    Zs2c(i)=y-mu2c(i)/(1-mu2c(i))*Z2c(i);
    Z2d(i)=y*log(1/mu2d(i));
    Zs2d(i)=y-mu2d(i)/(1-mu2d(i))*Z2d(i);
end

figure(1);
hold on;
plot(t,Z1a,'r');
plot(t,Zs1a,'g');
title('Problem 1a');
xlabel('Time (Gyr)');
ylabel('Metallicity');
legend('Gas','Stars');
plot(t,zeros(1,n),'--');
plot(t,ones(1,n)*0.1*Zsun);

figure(2);
hold on;
plot(t,Z1b,'r');
plot(t,Zs1b,'g');
title('Problem 1b');
xlabel('Time (Gyr)');
ylabel('Metallicity');
legend('Gas','Stars');
plot(t,zeros(1,n),'--');
plot(t,ones(1,n)*0.1*Zsun);

figure(3);
hold on;
plot(t,Z1c,'r');
plot(t,Zs1c,'g');
title('Problem 1c');
xlabel('Time (Gyr)');
ylabel('Metallicity');
legend('Gas','Stars');
plot(t,zeros(1,n),'--');
plot(t,ones(1,n)*0.1*Zsun);

figure(4);
hold on;
plot(t,Z1d,'r');
plot(t,Zs1d,'g');
title('Problem 1d');
xlabel('Time (Gyr)');
ylabel('Metallicity');
legend('Gas','Stars');
plot(t,zeros(1,n),'--');
plot(t,ones(1,n)*0.1*Zsun);

figure(5);
hold on;
plot(t,Z2a,'r');
plot(t,Zs2a,'g');
title('Problem 2a');
xlabel('Time (Gyr)');
ylabel('Metallicity');
legend('Gas','Stars');
plot(t,zeros(1,n),'--');
plot(t,ones(1,n)*0.1*Zsun);

figure(6);
hold on;
plot(t,Z2b,'r');
plot(t,Zs2b,'g');
title('Problem 2b');
xlabel('Time (Gyr)');
ylabel('Metallicity');
legend('Gas','Stars');
plot(t,zeros(1,n),'--');
plot(t,ones(1,n)*0.1*Zsun);

figure(7);
hold on;
plot(t,Z2c,'r');
plot(t,Zs2c,'g');
title('Problem 2c');
xlabel('Time (Gyr)');
ylabel('Metallicity');
legend('Gas','Stars');
line([1 1],get(gca,'ylim'));
line([2 2],get(gca,'ylim'));
line([3 3],get(gca,'ylim'));
line([4 4],get(gca,'ylim'));
line([5 5],get(gca,'ylim'));
line([6 6],get(gca,'ylim'));
plot(t,zeros(1,n),'--');
plot(t,ones(1,n)*0.1*Zsun);

figure(8);
hold on;
plot(t,Z2d,'r');
plot(t,Zs2d,'g');
title('Problem 2d');
xlabel('Time (Gyr)');
ylabel('Metallicity');
legend('Gas','Stars');
line([3 3],get(gca,'ylim'));
plot(t,zeros(1,n),'--');
plot(t,ones(1,n)*0.1*Zsun);