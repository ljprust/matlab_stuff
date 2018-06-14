% Logan Prust - Aer E 311 - Homework 3 - Problem 4

% Part (a)

clear,clc;
Pt=6; % MPa
Ae=3; % exit area
Ath=3+2.8*(-6)/6*(2+(-6)/6); % throat area
epsilon=10^(-6); % sets convergence criterion
x_initial=[0.001,2]; % sets initial values to start search
gamma=1.4;

% Newton's method for Me
for j=1:2
    n=1; % initializes the iteration counter
    error=1; % sets some arbitrary initial value for the error
    x=x_initial(j); % sets x equal to one of its initial values
while error>=epsilon
    n=n+1; % increases the iteration counter
    f=1/x*(2/(gamma+1)*(1+(gamma-1)/2*x^2))^((gamma+1)/2/(gamma-1))-Ae/Ath; % evaluates f(x)
    fprime=-1/(x^2)*(2/(gamma+1)*(1+(gamma-1)/2*x^2))^((gamma+1)/2/(gamma-1))+(gamma+1)/2/(gamma-1)*1/x*(2/(gamma+1)*(1+(gamma-1)/2*x^2))^((gamma+1)/2/(gamma-1)-1)*2*(gamma-1)/(gamma+1)*x;
    % ^ evaluates f'(x)
    x2=x-f/fprime; % evaluates x(n+1)
    error=abs(x2-x); % evaluates the absolute error
    x=x2; % makes x(n+1) the new x(n)
end
xvector(j)=x; % puts the two M values found into vector form
end
Me=min(xvector);

Patm1=Pt*(1+(gamma-1)/2*Me^2)^(-gamma/(gamma-1));

% Part (b)

Me=max(xvector);
Patm3=Pt*(1+(gamma-1)/2*Me^2)^(-gamma/(gamma-1));

% Part (c)

M1=Me;
M2=sqrt((M1^2+2/(gamma-1))/(2*gamma/(gamma-1)*M1^2-1));
Pratio=(1+gamma*M1^2)/(1+gamma*M2^2);
PtoverP1=(1+(gamma-1)/2*Me^2)^(gamma/(gamma-1));
Patm2=Pt*Pratio/PtoverP1;

% Part (d)

Patm=(Pt+Patm1)/2;
Me=sqrt(2/(gamma-1)*((Pt/Patm)^((gamma-1)/gamma)-1));
Astar=Ae*Me*(2/(gamma+1)*(1+(gamma-1)/2*Me^2))^((-gamma-1)/2/(gamma-1));

% Newton's method for throat Mach number
for j=1:2
    n=1; % initializes the iteration counter
    error=1; % sets some arbitrary initial value for the error
    x=x_initial(j); % sets x equal to one of its initial values
while error>=epsilon
    n=n+1; % increases the iteration counter
    f=1/x*(2/(gamma+1)*(1+(gamma-1)/2*x^2))^((gamma+1)/2/(gamma-1))-Ath/Astar; % evaluates f(x)
    fprime=-1/(x^2)*(2/(gamma+1)*(1+(gamma-1)/2*x^2))^((gamma+1)/2/(gamma-1))+(gamma+1)/2/(gamma-1)*1/x*(2/(gamma+1)*(1+(gamma-1)/2*x^2))^((gamma+1)/2/(gamma-1)-1)*2*(gamma-1)/(gamma+1)*x;
    % ^ evaluates f'(x)
    x2=x-f/fprime; % evaluates x(n+1)
    error=abs(x2-x); % evaluates the absolute error
    x=x2; % makes x(n+1) the new x(n)
end
xvector(j)=x; % puts the two M values found into vector form
end
Mth=min(xvector);

Pth=Pt*(1+(gamma-1)/2*Mth^2)^(-gamma/(gamma-1));

% Part (e)

Patm=(Patm1+Patm2)/2;
Me=sqrt(1/(gamma-1)*(-1+sqrt(1+2*(gamma-1)*(2/(gamma+1))^((gamma+1)/(gamma-1))*(Pt*Ath/Patm/Ae)^2)));
Aestar=Ae*Me*(2/(gamma+1)*(1+(gamma-1)/2*Me^2))^((-gamma-1)/2/(gamma-1));
Pte=Pt*Ath/Aestar;

% Newton's method for M1
for j=1:2
    n=1; % initializes the iteration counter
    error=1; % sets some arbitrary initial value for the error
    x=x_initial(j); % sets x equal to one of its initial values
while error>=epsilon
    n=n+1; % increases the iteration counter
    f=((gamma+1)*x^2/((gamma-1)*x^2+2))^(gamma/(gamma-1))*((gamma+1)/(2*gamma*x^2-gamma+1))^(1/(gamma-1))-Pte/Pt;
    fprime=gamma/(gamma-1)*((gamma+1)*x^2/((gamma-1)*x^2+2))^(gamma/(gamma-1)-1)*(2*(gamma+1)*x/((gamma-1)*x^2+2)-2*(gamma-1)*(gamma+1)*x^3/((gamma-1)*x^2+2)^2)*((gamma+1)/(2*gamma*x^2-gamma+1))^(1/(gamma-1))+((gamma+1)*x^2/((gamma-1)*x^2+2))^(gamma/(gamma-1))*1/(gamma-1)*((gamma+1)/(2*gamma*x^2-gamma+1))^(1/(gamma-1)-1)*(-(gamma+1)/(2*gamma*x^2-gamma+1)^2*4*gamma*x);
    x2=x-f/fprime; % evaluates x(n+1)
    error=abs(x2-x); % evaluates the absolute error
    x=x2; % makes x(n+1) the new x(n)
end
xvector(j)=x; % puts the two M values found into vector form
end
M1=xvector(2);

Ashock=Ath/M1*(2/(gamma+1)*(1+(gamma-1)/2*M1^2))^((gamma+1)/2/(gamma-1));
xshock=6*(-1+sqrt(1-0.357*(3-Ashock)));

% Part (f)

Patm=(Patm2+Patm3)/2;

% Newton's method for exit Mach number
for j=1:2
    n=1; % initializes the iteration counter
    error=1; % sets some arbitrary initial value for the error
    x=x_initial(j); % sets x equal to one of its initial values
while error>=epsilon
    n=n+1; % increases the iteration counter
    f=1/x*(2/(gamma+1)*(1+(gamma-1)/2*x^2))^((gamma+1)/2/(gamma-1))-Ae/Ath; % evaluates f(x)
    fprime=-1/(x^2)*(2/(gamma+1)*(1+(gamma-1)/2*x^2))^((gamma+1)/2/(gamma-1))+(gamma+1)/2/(gamma-1)*1/x*(2/(gamma+1)*(1+(gamma-1)/2*x^2))^((gamma+1)/2/(gamma-1)-1)*2*(gamma-1)/(gamma+1)*x;
    % ^ evaluates f'(x)
    x2=x-f/fprime; % evaluates x(n+1)
    error=abs(x2-x); % evaluates the absolute error
    x=x2; % makes x(n+1) the new x(n)
end
xvector(j)=x; % puts the two M values found into vector form
end
Me=max(xvector);

Pe=Pt*(1+(gamma-1)/2*Me^2)^(-gamma/(gamma-1));
sigma=asind(sqrt(1/2/gamma/Me^2*(gamma-1+(gamma+1)*Patm/Pe)));
delta=atand(2*cotd(sigma)*(Me^2*sind(sigma)^2-1)/(Me^2*(gamma+cosd(2*sigma))+2));
M2=1/sind(sigma-delta)*sqrt((Me^2*sind(sigma)^2+2/(gamma-1))/(2*gamma/(gamma-1)*Me^2*sind(sigma)^2-1));

% Part (g)

Patm=0.5*Patm3;
M1=Me;
theta1=asind(1/M1);
C=-theta1*3.14159/180-sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)/(gamma+1)*(M1^2-1)));
M2=sqrt(2/(gamma-1)*((Pt/Patm)^((gamma-1)/gamma)-1));
delta=(asin(1/M2)+C+sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)/(gamma+1)*(M2^2-1))))*180/3.14159;
theta2=asind(1/M2)-delta;