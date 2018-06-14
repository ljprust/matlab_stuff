% Logan Prust - Aer E 311 - Homework 5 - Problem 3

clear,clc;

% set constants
M=2.13;
x=linspace(0,1,1000);
alpha=[0,4,8];
alpharad=3.14159265358979323846/180*alpha;

% calculate pressure coefficients
for i=1:3
    for j=1:1000
        Cpplus(i,j)=2/sqrt(M^2-1)*(-alpharad(i)+0.2*(1-2*x(j)));
        Cpminus(i,j)=2/sqrt(M^2-1)*(alpharad(i)+0.2*(1-2*x(j)));
    end
end

% create plots
figure(1);
plot(x,Cpplus(1,:),'--');
hold on;
plot(x,Cpminus(1,:),'-');
xlabel('x*/c');
ylabel('Cp');
title('Cp versus x at alpha = 0 degrees');
legend('Upper Surface','Lower Surface');

figure(2);
plot(x,Cpplus(2,:),'--');
hold on;
plot(x,Cpminus(2,:),'-');
xlabel('x*/c');
ylabel('Cp');
title('Cp versus x at alpha = 4 degrees');
legend('Upper Surface','Lower Surface');

figure(3);
plot(x,Cpplus(3,:),'--');
hold on;
plot(x,Cpminus(3,:),'-');
xlabel('x*/c');
ylabel('Cp');
title('Cp versus x at alpha = 8 degrees');
legend('Upper Surface','Lower Surface');