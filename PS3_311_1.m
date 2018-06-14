% Logan Prust - Aer E 311 - Homework 3 - Problem 1

clear;

% get input parameters
M1=input('Upstream Mach Number: ');
deltascalar=input('Turning Angle (degrees): ');
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

sigmafinal=min(sigmavector); % gets the weak shock solution
fprintf('%20s %.15f %s \n','sigma =',sigmafinal,'degrees');

% calculate output parameters
M1n=M1*sind(sigmafinal);
M2n=sqrt((M1n^2+2/(gamma-1))/(2*gamma/(gamma-1)*M1n^2-1));
M2=M2n/sind(sigmafinal-deltascalar);
Ptratio=(1+gamma*M1n^2)/(1+gamma*M2n^2)*((1+(gamma-1)/2*M2n^2)/(1+(gamma-1)/2*M1n^2))^(gamma/(gamma-1));

% print output parameters
fprintf('\n %20s %.15f %s \n','M2 =',M2,'');
fprintf('%20s %.15f %s \n','Pt ratio =',Ptratio,'Pa');