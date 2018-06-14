% Logan Prust - Aer E 311 - Homework 5 - Problem 2

clear,clc;

M1=2;
P1=26510;
gamma=1.4;
deltascalar=10;
Pt=P1*(1+(gamma-1)/2*M1^2)^(gamma/(gamma-1));
n=100;
M2=linspace(1,10,n);
deltagrid=zeros(1,n);
theta1=asind(1/M1);
C=-theta1-sqrt((gamma+1)/(gamma-1))*atand(sqrt((gamma-1)/(gamma+1)*(M1^2-1)));

% find a value of delta for each value of M2
for i=1:n
    deltagrid(i)=C+asind(1/M2(i))+sqrt((gamma+1)/(gamma-1))*atand(sqrt((gamma-1)/(gamma+1)*(M2(i)^2-1)));
end

j=1;
while deltascalar>deltagrid(j)
    j=j+1;
end
M2initial=(deltascalar-deltagrid(j-1))/(deltagrid(j)-deltagrid(j-1))*(M2(j)-M2(j-1))+M2(j-1);

% use Newton's method to find M2 for each delta
epsilon=10^(-6);
deltanewton=linspace(0,20,n);
for k=1:n
    l=1; % initializes the iteration counter
    error=1; % sets some arbitrary initial value for the error
    M2scalar=M2initial; % sets M2 equal to its initial value
while error>=epsilon
    l=l+1; % increases the iteration counter
    f=sqrt(M2scalar^2-1)+sqrt((gamma+1)/(gamma-1))*tand(sqrt((gamma-1)/(gamma+1))*(asind(1/M2scalar)-deltanewton(k)+C));
    fprime=1/sqrt(M2scalar^2-1)*(M2scalar-1/M2scalar*(secd(sqrt((gamma-1)/(gamma+1))*(asind(1/M2scalar)-deltanewton(k)+C)))^2);
    M2scalar2=M2scalar-f/fprime; % evaluates M2(n+1)
    error=abs(M2scalar2-M2scalar); % evaluates the absolute error
    M2scalar=M2scalar2; % makes M2(n+1) the new M2(n)
end
    % calculate pressure ratio at each delta
    Pratio(k)=Pt/P1/(1+(gamma-1)/2*M2scalar^2)^(gamma/(gamma-1));
    M2vector(k)=M2scalar; % puts the M2 values found into vector form
end
P2upper=Pratio/600*P1;

% get input parameters
fprintf('\n');

% set some initial parameters
mu=asind(1/M1);
epsilon=10^(-6); % sets convergence criterion
sigma_initial=[mu+1,89]; % sets initial values to start search
alpha=linspace(0,20,100);

for k=1:100

for j=1:2
    n=1; % initializes the iteration counter
    error=1; % sets some arbitrary initial value for the error
    sigma=sigma_initial(j); % sets sigma equal to one of its initial values
while error>=epsilon
    n=n+1; % increases the iteration counter
    f=2*cotd(sigma)*(M1^2*(sind(sigma))^2-1)/(M1^2*(gamma+cosd(2*sigma))+2)-tand(alpha(k)); % evaluates f(sigma)
    fprime=-2*(cscd(sigma))^2*(M1^2*(sind(sigma))^2-1)/(M1^2*(gamma+cosd(2*sigma))+2)+2*cotd(sigma)*(2*M1^2*sind(sigma)*cosd(sigma))/(M1^2*(gamma+cosd(2*sigma))+2)-2*cotd(sigma)*(M1^2*(sind(sigma))^2-1)/(M1^2*(gamma+cosd(2*sigma))+2)^2*(-2*M1^2*sind(2*sigma));
    % ^ evaluates f'(sigma)
    sigma2=sigma-f/fprime; % evaluates sigma(n+1)
    error=abs(sigma2-sigma); % evaluates the absolute error
    sigma=sigma2; % makes sigma(n+1) the new sigma(n)
end
sigmavector(j)=sigma; % puts the two sigma values found into vector form
end

% chooses the right sigma value for the strong or weak condition
strongorweak=' weak ';
if strongorweak=='strong'
    sigmafinal=max(sigmavector);
else
    sigmafinal=min(sigmavector);
end

% calculate output parameters
M1n=M1*sind(sigmafinal);
M2n=sqrt((M1n^2+2/(gamma-1))/(2*gamma/(gamma-1)*M1n^2-1));
M2=M2n/sind(sigmafinal-deltascalar);
P2lower(k)=P1/600*(1+gamma*M1n^2)/(1+gamma*M2n^2);

end

rhoinf=0.4135;
Tinf=223.3;
ainf=sqrt(gamma*287*Tinf);
vinf=M1*ainf;
for z=1:100
    Cl(z)=(P2lower(z)-P2upper(z))*2/rhoinf/vinf*cosd(alpha(z));
    Cd(z)=(P2lower(z)-P2upper(z))*2/rhoinf/vinf*sind(alpha(z));
    Clthin(z)=4*alpha(z)*3.14159265358979323846/180/sqrt(3);
    Cdthin(z)=4*(alpha(z)*3.14159265358979323846/180)^2/sqrt(3);
end

figure(1);
plot(alpha,Cl,'-');
hold on;
plot(alpha,Clthin,'--');
xlabel('Angle of Attack (degrees)');
ylabel('Coefficients of Lift');
title('Lift Coefficients vs Angle of Attack');
legend('Shock/Expansion','Thin Airfoil');

figure(2);
plot(alpha,Cd,'-');
hold on;
plot(alpha,Cdthin,'--');
xlabel('Angle of Attack (degrees)');
ylabel('Coefficients of Drag');
title('Drag Coefficients vs Angle of Attack');
legend('Shock/Expansion','Thin Airfoil');