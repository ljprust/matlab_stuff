% Logan Prust - Astro 586 - Homework 5 - Problem 3(c)

clear,clc;

e=0.2;
t=linspace(0,4*pi,100);

for i=1:100
    xkep(i)=e*cos(t(i));
    ykep(i)=-2*e*sin(t(i));
    xflat(i)=e*cos(sqrt(2)*t(i));
    yflat(i)=-sqrt(2)*e*sin(sqrt(2)*t(i));
    rkep(i)=1+xkep(i);
    thetakep(i)=t(i)+ykep(i);
    rflat(i)=1+xflat(i);
    thetaflat(i)=t(i)+yflat(i);
end

polar(thetakep,rkep);
hold on
polar(thetaflat,rflat,'r');
title('Comparison of Epicycles');
legend('Keplerian','Logarithmic');