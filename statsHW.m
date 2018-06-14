clear,clc;
decays=[0,linspace(1,1,9),linspace(2,2,20),linspace(3,3,24),linspace(4,4,19),linspace(5,5,11),linspace(6,6,11),8,8,8,9,10];
hist(decays,0:10);
hold on

muparent=3.7;
muobserved=3.61;
parent0=exp(-1*muparent)*100;
observed0=exp(-1*muobserved)*100;
plot(0,parent0,'r*');
plot(0,observed0,'g*');
for x=1:10
    parent(x)=muparent^(x)/factorial(x)*exp(-1*muparent)*100;
    observed(x)=muobserved^(x)/factorial(x)*exp(-1*muobserved)*100;
    plot(x,parent(x),'r*');
    plot(x,observed(x),'g*');
end
title('Decays Per Minute');