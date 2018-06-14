% Logan Prust - Aer E 311 - Homework 4 - Problem 1

clear;

% generates Mach numbers for plot
gamma=1.4;
M=linspace(0,4,1000);
% generates fR values for plot
for i=1:1000
    fR(i)=(1+(gamma-1)/2*M(i)^2)/(1+gamma*M(i)^2)^2*M(i)^2;
end
% creates plot
plot(M,fR)
xlabel('Mach Number');
ylabel('fR');
title('fR vs Mach Number');

clear;

% gets desired upstream Mach number
M1=input('M1 = ');
% sets some parameters for the problem
T1=1000;
h=4.4e7;
gamma=1.4;
R=287;
Cp=gamma*R/(gamma-1);
mfdot=8;
mdot=800;
% calculates relevant values
f=mfdot/mdot;
Tt1=T1*(1+(gamma-1)/2*M1^2);
Ttratio=1+f*h/Cp/Tt1;
fR=(1+(gamma-1)/2*M1^2)/(1+gamma*M1^2)^2*M1^2;
f2=Ttratio*fR;
% finds downstream Mach number for either subsonic or supersonic condition
if M1<1
    M2=sqrt((1-2*gamma*f2-sqrt((2*gamma*f2-1)^2-4*f2*(f2*gamma^2-(gamma-1)/2)))/(2*f2*gamma^2+1-gamma));
else
    M2=sqrt((1-2*gamma*f2+sqrt((2*gamma*f2-1)^2-4*f2*(f2*gamma^2-(gamma-1)/2)))/(2*f2*gamma^2+1-gamma));
end
Ptratio=((1+(gamma-1)/2*M2^2)/(1+(gamma-1)/2*M1^2))^(gamma/(gamma-1))*(1+gamma*M1^2)/(1+gamma*M2^2);
% prints out desired output values
fprintf('%20s %.15f %s \n','M2 =',M2,'');
fprintf('%20s %.15f %s \n','Pt2/Pt1 =',Ptratio,'');