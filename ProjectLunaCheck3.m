clear,clc;

digits=1000;
format long;

% % v0 in DU/TU, lambda1 in degrees
% v0=1.385
% lambda1=44.527355

% v0 = 1.39
% lambda1 = 0.541465362*180/pi

v0 = 1.372
lambda1 = 35.312

DU=6378.1; % (km/DU) conversion between canonical and SI
muESI=3.986e5; % (SI) mu of Earth
TU=sqrt(DU^3/muESI); % (sec/TU) conversion between canonical and SI
phi0=0; % initial flight path angle
RE=1; % (canonical) radius of Earth
alt0=185/DU; % (canonical) altitude of parking orbit
R0=RE+alt0;% (canonical) radius of parking orbit
% R0=1.05;
V0esc=sqrt(2/R0);% (canonical) escape speed of parking orbit
RsoiMSI=66183; % (SI) radius of lunar sphere of influence
RsoiM=RsoiMSI/DU; % (canonical) radius of lunar sphere of influence
DSI=384400; % (SI) Earth-Moon distance
D=DSI/DU; % (canonical) Earth-Moon distance
muE=1; % (canonical) mu of Earth
muM=4902.8; % (SI) mu of Moon
RM=1738; % (SI) radius of Moon
wMSI=2.64911e-6; % (rad/sec) orbital angular velocity of Moon
wM=wMSI*TU*180/pi; % (degrees/TU) orbital angular velocity of Moon

HT=R0*v0; % (canonical) same parameters as before, but to check the candidate trajectory
ET=v0^2/2-muE/R0; % (canonical)
pT=HT^2/muE; % (canonical)
aT=-muE/(2*ET); % (canonical)
eT=sqrt(1-pT/aT); % (dimensionless)
R1=sqrt(RsoiM^2+D^2-2*RsoiM*D*cosd(lambda1)); % (canonical)
nu1=acosd((pT-R1)/(R1*eT)); % (degrees)

hold on;
circle(0,0,1);
axis([-2 71 -28 28]);
circle(D,0,RM/DU);
circle(D,0,RsoiM);

if imag(nu1) ~= 0
    angmax = 2*pi;
    ang_t = 0:0.01:angmax;
    r_t = pT./(1+eT*cos(ang_t));
    phase_t = pi;
    ang_t = ang_t - phase_t;
    polar(ang_t,r_t,'r');
else
    angmax = nu1;
end

gamma1=asind(RsoiM/R1*sind(lambda1)); % (degrees)
% nu0=acosd((pT-R0)/(R0*eT)); % (degrees)
nu0=0;
EA0=acos((eT+cosd(nu0))/(1+eT*cosd(nu0))); % (radians)
EA1=acos((eT+cosd(nu1))/(1+eT*cosd(nu1))); % (radians)
TOF01=sqrt(aT^3/muE)*(EA1-eT*sin(EA1)-EA0+eT*sin(EA0)); % (canonical)
TOF01hr=TOF01*TU/3600; % (hours)
beta01=wM*TOF01; % (degrees)
gamma0=nu1-nu0-gamma1-wM*TOF01; % (degrees)

v1=sqrt(2*(ET+muE/R1)); % (canonical)
phi1=acosd(HT/R1/v1); % (degrees)

R2=RsoiMSI; % (SI)
vM=wMSI*DSI; % (SI)
% vM=sqrt(muESI/DSI);
v1SI=v1*DU/TU; % (SI)
v2=sqrt(v1SI^2+vM^2-2*v1SI*vM*cosd(phi1-gamma1)); % (SI)
eps2=asind(vM/v2*cosd(lambda1)-v1SI/v2*cosd(lambda1+gamma1-phi1)); % (degrees)
phi2=abs(eps2)-90;% (degrees)
EM=v2^2/2-muM/R2; % (SI)
HM=R2*v2*abs(sind(eps2)); % (SI)
pM=HM^2/muM; % (SI)
eM=sqrt(1+2*EM*HM^2/muM^2); % (dimensionless)
aM=pM/(1-eM^2); % (SI)
RpM=pM/(1+eM); % (SI)
nu2=acosd((pM-R2)/(R2*eM)); % (degrees)
% nu3=acosd((pM-RpM)/(RpM*eM)); % (degrees)
nu3=0;
EA2=acos((eM+cosd(nu2))/(1+eM*cosd(nu2))); % (radians)
EA3=acos((eM+cosd(nu3))/(1+eM*cosd(nu3))); % (radians)
TOF23=-sqrt(aM^3/muM)*(EA3-eM*sin(EA3)-EA2+eM*sin(EA2)); % (SI)
TOF23hr=TOF23/3600; % (hours)
TOFhr=TOF01hr+TOF23hr; % (hours)
altM=RpM-RM; % (SI)
vpM=sqrt(2*(EM+muM/RpM)); % (SI)

ang_t = 0:0.01:angmax;
r_t = pT./(1+eT*cosd(ang_t));
phase_t = beta01 + gamma0;
ang_t = ang_t - phase_t;
ang_t = ang_t*pi/180;
polar(ang_t,r_t,'r');

fprintf('%20s %.15f %s \n','Altitude =',altM,'km')
if altM>=119.9 && altM<=120.1
    disp('Check');
    disp(' ');
else
    disp('No Check');
    disp(' ');
end
fprintf('%20s %.15f %s \n','Energy =',ET,'DU^2/TU^2')
if ET<0
    disp('Check - Elliptical Cislunar Approach');
    disp(' ');
else
    disp('No Check');
    disp(' ');
end
fprintf('%20s %.15f %s \n','eps2 =',eps2,'degrees')
if eps2>0
    disp('No Check - Direct');
    disp(' ');
else
    disp('Check - Retrograde');
    disp(' ');
end
% Results in proper order
fprintf('Departure velocity: %35f DU \n',v0)
fprintf('Departure flight path angle: %26f deg \n',phi0)
fprintf('Lambda 1: %46f deg \n',lambda1)
fprintf('Energy at departure: %35f DU^2/TU^2 \n',ET)
fprintf('True anomaly angle at departure: %22f deg \n',nu0)
fprintf('True anomaly angle at SOI: %30f deg \n',nu1)
fprintf('Beta01: %48f deg \n',beta01)
fprintf('Angle gamma0: %43f deg \n\n',gamma0)
fprintf('Radius from Earth at SOI: %30f DU \n',R1)
fprintf('Velocity at SOI(Earth Frame): %25f DU/TU \n',v1)
fprintf('Flight path angle at SOI(Earth Frame): %17f deg \n',phi1)
fprintf('Angle gamma1: %41f deg \n',gamma1)
fprintf('Semi major axis of transfer ellipse: %19f DU \n',aT)
fprintf('Eccentricity of transfer ellipse: %21f\n\n',eT)
fprintf('Velocity at SOI(Lunar Frame): %25f km/s \n',v2)
% fprintf('Flight path angle at SOI(Lunar Frame): %17f deg \n',phi2)
fprintf('Angle epsilon2: %39f deg \n',eps2)
fprintf('Energy inside SOI: %36f km^2/s^2 \n',EM)
fprintf('Semi major axis of lunar orbit: %27f km \n',aM)
fprintf('Eccentricity of lunar orbit: %26f\n\n',eM)
fprintf('Radius of Periselenium: %34f km \n',RpM)
fprintf('Altitude at Periselenium: %32f km \n',altM)
fprintf('Velocity at Periselenium: %29f km/s \n\n',vpM)
fprintf('TOF from translunar injection to the Moon SOI: %9f hours \n',TOF01hr)
fprintf('TOF from SOI entry to periselenium: %20f hours \n',TOF23hr)
fprintf('Total TOF: %45f hours \n',TOFhr)

% plot trajectory

ang_M = 0:0.01:nu2;
r_M = pM/DU./(1+eM*cosd(ang_M));
phase_M = 180 - lambda1 - nu2;
ang_M = ang_M + phase_M;

if eps2>0
    ang_M = 2*(180 - lambda1) - ang_M;
end

x_M = r_M.*cosd(ang_M) + D;
y_M = r_M.*sind(ang_M);
plot(x_M,y_M,'r');