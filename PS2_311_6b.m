% Logan Prust - Aer E 311 - Homework 2 - Problem 6 - Part (b)

clear;

% get input parameters
M1=input('Upstream Mach Number: ');
P1=input('Upstream Pressure: ');
rho1=input('Upstream Density: ');
T1=input('Upstream Temperature: ');
Pt1=input('Upstream Total Pressure: ');
deltascalar=input('Turning Angle (degrees): ');
strongorweak=input('Strong or weak shock? ');
fprintf('\n');

% set some initial parameters
gamma=1.4;
mu=asind(1/M1);
epsilon=10^(-6); % sets convergence criterion
sigma_initial=[mu+1,89]; % sets initial values to start search

for j=1:2
    n=1; % initializes the iteration counter
    error=1; % sets some arbitrary initial value for the error
    sigma=sigma_initial(j); % sets sigma equal to one of its initial values
while error>=epsilon
    n=n+1; % increases the iteration counter
    f=2*cotd(sigma)*(M1^2*(sind(sigma))^2-1)/(M1^2*(gamma+cosd(2*sigma))+2)-tand(deltascalar); % evaluates f(sigma)
    fprime=-2*(cscd(sigma))^2*(M1^2*(sind(sigma))^2-1)/(M1^2*(gamma+cosd(2*sigma))+2)+2*cotd(sigma)*(2*M1^2*sind(sigma)*cosd(sigma))/(M1^2*(gamma+cosd(2*sigma))+2)-2*cotd(sigma)*(M1^2*(sind(sigma))^2-1)/(M1^2*(gamma+cosd(2*sigma))+2)^2*(-2*M1^2*sind(2*sigma));
    % ^ evaluates f'(sigma)
    sigma2=sigma-f/fprime; % evaluates sigma(n+1)
    error=abs(sigma2-sigma); % evaluates the absolute error
    sigma=sigma2; % makes sigma(n+1) the new sigma(n)
end
sigmavector(j)=sigma; % puts the two sigma values found into vector form
end

% chooses the right sigma value for the strong or weak condition
if strongorweak=='strong'
    sigmafinal=max(sigmavector);
else
    sigmafinal=min(sigmavector);
end
fprintf('%20s %.15f %s \n','sigma =',sigmafinal,'degrees');

% calculate output parameters
M1n=M1*sind(sigmafinal);
M2n=sqrt((M1n^2+2/(gamma-1))/(2*gamma/(gamma-1)*M1n^2-1));
M2=M2n/sind(sigmafinal-deltascalar);
P2=P1*(1+gamma*M1n^2)/(1+gamma*M2n^2);
T2=T1*(1+(gamma-1)/2*M1n^2)/(1+(gamma-1)/2*M2n^2);
rho2=rho1*M1n/M2n*sqrt((1+(gamma-1)/2*M2n^2)/(1+(gamma-1)/2*M1n^2));
Pt2=Pt1*(1+gamma*M1n^2)/(1+gamma*M2n^2)*((1+(gamma-1)/2*M2n^2)/(1+(gamma-1)/2*M1n^2))^(gamma/(gamma-1));

% print output parameters
fprintf('\n %20s %.15f %s \n','M2 =',M2,'');
fprintf('%20s %.15f %s \n','P2 =',P2,'Pa');
fprintf('%20s %.15f %s \n','T2 =',T2,'K');
fprintf('%20s %.15f %s \n','rho2 =',rho2,'kg/m^3');
fprintf('%20s %.15f %s \n','Pt2 =',Pt2,'Pa');

% set some parameters for the plots
n=1000;
M1=[1.5,2];
mu=[asind(1/1.5),asind(1/2)];

% create plots
for i=1:2
    sigma=linspace(mu(i),90,n);
    for j=1:n
        delta(j)=atand(2*cotd(sigma(j))*(M1(i)^2*(sind(sigma(j)))^2-1)/(M1(i)^2*(gamma+cosd(2*sigma(j)))+2));
    end
    if i==1
        plot(delta,sigma,'b');
    else
        plot(delta,sigma,'g');
    end
    hold on;
end

% make plots look pretty
plot([deltascalar,deltascalar],sigmavector,'r*');
legend('M = 1.5','M = 2','Solution Points');
axis([0 25 20 100]);
title('Shock Angle vs Turning Angle');
xlabel('Turning Angle (degrees)');
ylabel('Shock Angle (degrees)');