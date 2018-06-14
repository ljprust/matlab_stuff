% Aer E 355 HW 7 Problem 2

clear,clc;

% M=0.2

W=5.640E5;
g=32.2;
m=W/g;
Xu=-3.661E2;
Xw=2.137E3;
theta0=0;
Zu=-3.538E3;
Zwdot=5.851E2;
Zw=-8.969E3;
Zq=-1.090E5;
u0=221;
Iy=3.23E7;
Mu=3.779E3;
Mwdot=-7.946E3;
Mw=-5.717E4;
Mq=-1.153E7;

Along=[Xu/m Xw/m 0 -g*cos(theta0);Zu/(m-Zwdot) Zw/(m-Zwdot) (Zq+m*u0)/(m-Zwdot) -m*g*sin(theta0)/(m-Zwdot);(Mu+Mwdot*Zu/(m-Zwdot))/Iy (Mw+Mwdot*Zw/(m-Zwdot))/Iy (Mq+Mwdot*(Zq+m*u0)/(m-Zwdot))/Iy -Mwdot*m*g*sin(theta0)/Iy/(m-Zwdot);0 0 1 0];
case1eig=eig(Along)

% M=0.5

W=6.366E5;
m=W/g;
Xu=-4.883E1;
Xw=1.546E3;
Zu=-1.342E3;
Zwdot=3.104E2;
Zw=-8.561E3;
Zq=-1.263E5;
u0=518;
Iy=3.31E7;
Mu=8.176E3;
Mwdot=-4.138E3;
Mw=-5.627E4;
Mq=-1.394E7;

Along=[Xu/m Xw/m 0 -g*cos(theta0);Zu/(m-Zwdot) Zw/(m-Zwdot) (Zq+m*u0)/(m-Zwdot) -m*g*sin(theta0)/(m-Zwdot);(Mu+Mwdot*Zu/(m-Zwdot))/Iy (Mw+Mwdot*Zw/(m-Zwdot))/Iy (Mq+Mwdot*(Zq+m*u0)/(m-Zwdot))/Iy -Mwdot*m*g*sin(theta0)/Iy/(m-Zwdot);0 0 1 0];
case2eig=eig(Along)

% M=0.9

W=6.366E5;
m=W/g;
Xu=-3.954E2;
Xw=3.144E2;
Zu=-8.383E2;
Zwdot=1.214E2;
Zw=-7.928E3;
Zq=-1.327E5;
u0=871;
Iy=3.31E7;
Mu=-2.062E3;
Mwdot=-5.296E3;
Mw=-6.289E4;
Mq=-1.327E7;

Along=[Xu/m Xw/m 0 -g*cos(theta0);Zu/(m-Zwdot) Zw/(m-Zwdot) (Zq+m*u0)/(m-Zwdot) -m*g*sin(theta0)/(m-Zwdot);(Mu+Mwdot*Zu/(m-Zwdot))/Iy (Mw+Mwdot*Zw/(m-Zwdot))/Iy (Mq+Mwdot*(Zq+m*u0)/(m-Zwdot))/Iy -Mwdot*m*g*sin(theta0)/Iy/(m-Zwdot);0 0 1 0];
case3eig=eig(Along)