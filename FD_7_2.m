% Aer E 355 HW 7 Problem 5

clear,clc;

CL=.192;
W=10675;
g=9.81;
m=W/g;
S=14.9;
b=9.14;
u0=77.3;
roh=1.25;
theta0=0;

Ix=230;
Iz=1778;
Izx=0;

Ixprime=(Ix*Iz-Izx^2)/Iz;
Izprime=(Ix*Iz-Izx^2)/Ix;
Izxprime=Izx/(Ix*Iz-Izx^2);

Cybeta=-.14;
Cyp=-.039;
Cyr=.165;

Clbeta=-.0689*CL-.0917*CL;
Clp=-.441;
Clr=-.0144+.271*CL;

Cnbeta=.01326+.017*CL^2;
Cnp=-.00109-.0966*CL;
Cnr=-.048-.0238*CL^2;

Yv=.5*roh*u0*S*Cybeta;
Yp=.25*roh*u0*b*S*Cyp;
Yr=.25*roh*u0*b*S*Cyr;

Lv=.5*roh*u0*b*S*Clbeta;
Lp=.25*roh*u0*b^2*S*Clp;
Lr=.25*roh*u0*b^2*S*Clr;

Nv=.5*roh*u0*b*S*Cnbeta;
Np=.25*roh*u0*b^2*S*Cnp;
Nr=.25*roh*u0*b^2*S*Cnr;

Alat=[Yv/m Yp/m Yr/m-u0 g*cos(theta0);Lv/Ixprime+Izxprime*Nv Lp/Ixprime+Izxprime*Np Lr/Ixprime+Izxprime*Nr 0;Izxprime*Lv+Nv/Izprime Izxprime*Lp+Np/Izprime Izxprime*Lr+Nr/Izprime 0;0 1 tan(theta0) 0]

lambda=eig(Alat)