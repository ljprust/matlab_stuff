% Logan Prust - Aer E 311 - Homework 2 - Problem 7

clear;

% get input parameters
M1=input('Upstream Mach Number: ');
P1=input('Upstream Pressure: ');
deltascalar=input('Turning Angle (degrees): ');

% set some initial parameters for the code
gamma=1.4;
Pt=P1*(1+(gamma-1)/2*M1^2)^(gamma/(gamma-1));
n=1000;
M2=linspace(1,10,n);
deltagrid=zeros(1,n);
theta1=asind(1/M1);
C=-theta1-sqrt((gamma+1)/(gamma-1))*atand(sqrt((gamma-1)/(gamma+1)*(M1^2-1)));

% find a value of delta for each value of M2
for i=1:n
    deltagrid(i)=C+asind(1/M2(i))+sqrt((gamma+1)/(gamma-1))*atand(sqrt((gamma-1)/(gamma+1)*(M2(i)^2-1)));
end

% use the grid of delta values to interpolate an inital value for M2
j=1;
while deltascalar>deltagrid(j)
    j=j+1;
end
M2initial=(deltascalar-deltagrid(j-1))/(deltagrid(j)-deltagrid(j-1))*(M2(j)-M2(j-1))+M2(j-1);

% use Newton's method to find M2 for each delta
epsilon=10^(-6);
deltanewton=linspace(0,40,n);
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

% find the values of M2 and P2 at the specified point
error=1;
M2point=M2initial;
while error>=epsilon
    f=sqrt(M2point^2-1)+sqrt((gamma+1)/(gamma-1))*tand(sqrt((gamma-1)/(gamma+1))*(asind(1/M2point)-deltascalar+C));
    fprime=1/sqrt(M2point^2-1)*(M2point-1/M2point*(secd(sqrt((gamma-1)/(gamma+1))*(asind(1/M2point)-deltascalar+C)))^2);
    M2point2=M2point-f/fprime;
    error=abs(M2point2-M2point);
    M2point=M2point2;
end
Ppoint=Pt/P1/(1+(gamma-1)/2*M2point^2)^(gamma/(gamma-1));
P2=Ppoint*P1;

% prints output parameters
fprintf('\n %20s %.15f %s \n','M2 =',M2point,'');
fprintf('%20s %.15f %s \n','P2 =',P2,'Pa');

% construct plots
plot(deltanewton,M2vector,'g');
hold on;
plot(deltanewton,Pratio,'b');
plot(deltascalar,M2point,'r*');
legend('M2','P2/P1','Solution Points');
plot(deltascalar,Ppoint,'r*');
title('Prandtl-Meyer Expansion Properties');
xlabel('Turning Angle (degrees)');
ylabel('M2 and P2/P1');