% Logan Prust - Aer E 311 - Homework 4 - Problem 2

clear,clc;

% Part (d) ------------------------------

% sets some initial parameters
gamma=1.4;
R=287;
epsilon=10^(-6);
mdot=970.04;
h0=10000;
kappa=1.5*10^(-4);
P0=100000;
Patm=P0*exp(-kappa*h0);
Tt=2500;
Pt=6000000;
Astar=0.2;
A3rd=3.15;

% uses Newton's method to find Me
Ae=linspace(A3rd/5,A3rd*5,1000);
for i=1:1000

Me_initial=[0.001,2];
for j=1:2
    n=1;
    error=1;
    Me=Me_initial(j);
while error>=epsilon
    n=n+1;
    f=1/Me*(2/(gamma+1)*(1+(gamma-1)/2*Me^2))^((gamma+1)/2/(gamma-1))-Ae(i)/Astar; % evaluates f(x)
    fprime=-1/(Me^2)*(2/(gamma+1)*(1+(gamma-1)/2*Me^2))^((gamma+1)/2/(gamma-1))+(gamma+1)/2/(gamma-1)*1/Me*(2/(gamma+1)*(1+(gamma-1)/2*Me^2))^((gamma+1)/2/(gamma-1)-1)*2*(gamma-1)/(gamma+1)*Me;
    Me2=Me-f/fprime;
    error=abs(Me2-Me);
    Me=Me2;
end
Mevector(j)=Me;
end
Mexit(i)=max(Mevector); % the supersonic solution is selected

% some parameters are calculated for the plots
Te=Tt/(1+(gamma-1)/2*Mexit(i)^2);
Pe=Pt*(Te/Tt)^(gamma/(gamma-1));
ue(i)=Mexit(i)*sqrt(gamma*R*Te);
F(i)=mdot*ue(i)+(Pe-Patm)*Ae(i);

end

% creates the plots
figure(1);
plot(Ae,F,A3rd,1.94e6,'r*');
title('Thrust vs Exit Area');
xlabel('Exit Area (m^2)');
ylabel('Thrust (N)');
figure(2);
plot(Ae,ue,A3rd,2002,'r*');
title('Exit Velocity vs Exit Area');
xlabel('Exit Area (m^2)');
ylabel('Exit Velocity (m/s)');

% Part (f) ------------------------------

% uses the 3rd critical conditions
Amax=A3rd;
h=linspace(0,50,1000);
ue3rd=2002;
Pe3rd=22300;

% calculates thrust at all altitudes
for k=1:1000
    F(k)=mdot*ue3rd+(Pe3rd-P0*exp(-kappa*h(k)*1000))*Amax;
end

% creates the plot
figure(3);
plot(F,h);
title('Thrust vs Altitude');
xlabel('Thrust (N)');
ylabel('Altitude (km)');