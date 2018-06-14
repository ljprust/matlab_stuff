% Logan Prust - Aer E 331 - Homework 6

% PROBLEM 3A

clear;

kvec=0:10;
u=zeros(1,length(kvec));
u(1)=1;
y=zeros(1,length(kvec));

y(1+1)=u(1+1)+0.5*u(1-1+1)+y(1-1+1)-0.09*0;

for k=2:10
    j=k+1;
    y(j)=u(j)+0.5*u(j-1)+y(j-1)-0.09*y(j-2);
end

figure(1);
plot(kvec,y);
xlabel('k (Discrete Time)');
ylabel('y');
title('System Impulse Response');

% PROBLEM 4B

clear;

kvec=0:10;
u=ones(1,length(kvec));
y=zeros(1,length(kvec));

y(1+1)=2*u(1+1-1)+3*0-2*0;

for k=2:10
    j=k+1;
    y(j)=2*u(j-1)+3*y(j-1)-2*y(j-2);
end

figure(2);
plot(kvec,y);
xlabel('k (Discrete Time)');
ylabel('y');
title('System Step Response');