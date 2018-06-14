% Logan Prust - Astro 586 - Exam 1

clear,clc;

n=100;
t=linspace(0,10,n);
Zsun=0.02;
y=Zsun;
taub=1;
taud=10;

for i=1:n
    mub(i)=0.99*exp(-t(i)/taub)+0.01;
    Zb(i)=y*log(1/mub(i));
    Zsb(i)=y*(1-mub(i)/(1-mub(i))*log(1/mub(i)));
    mud(i)=1/(1+0.9*exp(-t(i)/taud)*t(i));
    nud(i)=0.9*exp(-t(i)/taud)*t(i);
    Zd(i)=y*(1-exp(-nud(i)));
    Zsd(i)=y-mud(i)/(1-mud(i))*Zd(i);
end

figure(1);
hold on;
plot(t,Zb,'r');
plot(t,Zsb,'b');
title('Figure 3');
xlabel('Time (Gyr)');
ylabel('Metallicity');
legend('Gas (Newly Formed Stars)','Stars');

figure(2);
hold on;
plot(t,Zd,'r');
plot(t,Zsd,'b');
title('Figure 4');
xlabel('Time (Gyr)');
ylabel('Metallicity');
legend('Gas (Newly Formed Stars)','Stars');