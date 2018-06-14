% Logan Prust - Aer E 351 - Project Luna
clear,clc;
M=1000; % number of velocity iterations
N=1000; % number of lambda iterations
v0=linspace(1.382,1.383,M); % set up v0 array
lambda1=linspace(45,50,N); % set up lambda1 array
RpM=zeros(M,N); % preallocate matrices
altM=zeros(M,N);
vpM=zeros(M,N);
A=zeros(M,N);
B=zeros(M,N);
C=zeros(M,N);
DU=6378.1; % (km/DU) conversion between canonical and SI
muESI=3.9860044e5; % (SI) mu of Earth
TU=sqrt(DU^3/muESI); % (sec/TU) conversion between canonical and SI
phi0=0; % initial flight path angle
RE=1; % (canonical) radius of Earth
alt0=185/DU; % (canonical) altitude of parking orbit
R0=RE+alt0; % (canonical) radius of parking orbit
RsoiMSI=66183; % (SI) radius of lunar sphere of influence
RsoiM=66183/DU; % (canonical) radius of lunar sphere of influence
DSI=384400; % (SI) Earth-Moon distance
D=384400/DU; % (canonical) Earth-Moon distance
muE=1; % (canonical) mu of Earth
muM=4.9028e3; % (SI) mu of Moon
RM=1738; % (SI) radius of Moon
wMSI=2.6491065012e-6; % (rad/sec) orbital angular velocity of Moon
wM=(2.6491065012e-6)*TU*180/pi; % (degrees/TU) orbital angular velocity of Moon
for i=1:M % iterate over velocity
    for j=1:N % iterate over lambda
HT=R0*v0(i); % (canonical) angular momentum of elliptical orbit
ET=v0(i)^2/2-muE/R0; % (canonical) energy of elliptical orbit
pT=HT^2/muE; % (canonical) semi-latus parameter of elliptical orbit
aT=-muE/(2*ET); % (canonical) semi-major axis of elliptical orbit
eT=sqrt(1-pT/aT); % (dimensionless) eccentricity of elliptical orbit
R1=sqrt(RsoiM^2+D^2-2*RsoiM*D*cosd(lambda1(j))); % (canonical) radius from Earth at SOI
v1=sqrt(2*(ET+muE/R1)); % (canonical) velocity in Earth frame at SOI
phi1=acosd(real(HT/R1/v1)); % (degrees) flight-path angle in Earth frame at SOI
gamma1=asind(RsoiM/R1*sind(lambda1(j))); % (degrees) gamma at SOI
nu0=real(acosd((pT-R0)/(R0*eT))); % (degrees) inital true anomaly
nu1=real(acosd((pT-R1)/(R1*eT))); % (degrees) true anomaly in Earth frame at SOI
EA0=acos((eT+cosd(nu0))/(1+eT*cosd(nu0))); % (radians) initial eccentric anomaly
EA1=acos((eT+cosd(nu1))/(1+eT*cosd(nu1))); % (radians) eccentric anomaly in Earth frame at SOI
TOF01=sqrt(aT^3/muE)*(EA1-eT*sin(EA1)-EA0+eT*sin(EA0)); % (canonical) time of flight up to SOI
TOF01hr=TOF01*TU/3600; % (hours) time of flight up to SOI
gamma0=nu1-nu0-gamma1-wM*TOF01; % (degrees) lead angle
R2=RsoiMSI; % (SI) radius from Moon at SOI
vM=wMSI*DSI; % (SI) orbital velocity of the Moon
v1SI=v1*DU/TU; % (SI) velocity in Earth frame at SOI
v2=sqrt(v1SI^2+vM^2-2*v1SI*vM*cosd(real(phi1-gamma1))); % (SI) velocity in Moon frame at SOI
eps2=asind(real(vM/v2*cosd(lambda1(j))-v1SI/v2*cosd(real(lambda1(j)+gamma1-phi1)))); % (degrees) epsilon at SOI
EM=v2^2/2-muM/R2; % (SI) energy in Moon frame
HM=R2*v2*sind(real(eps2)); % (SI) angular momentum in Moon frame
pM=HM^2/muM; % (SI) semi-latus parameter in Moon frame
eM=sqrt(1+2*EM*HM^2/muM^2); % (dimensionless) eccentricity in Moon frame
aM=pM/(1-eM^2); % (SI) semi-major axis in Moon frame
RpM(i,j)=pM/(1+eM); % (SI) periselenium radius
nu2=real(acosd(real((pM-R2)/(R2*eM)))); % (degrees) true anomaly in Moon frame at SOI
nu3=real(acosd(real((pM-RpM(i,j))/(RpM(i,j)*eM)))); % (degrees) true anomaly in Moon frame at periselenium
EA2=acos((eM+cosd(nu2))/(1+eM*cosd(nu2))); % (radians) eccentric anomaly in Moon frame at SOI
EA3=acos((eM+cosd(nu3))/(1+eM*cosd(nu3))); % (radians) eccentric anomaly in Moon frame at periselenium
TOF23=-sqrt(aM^3/muM)*(EA3-eM*sin(EA3)-EA2+eM*sin(EA2)); % (SI) time of flight from SOI to periselenium
TOF23hr=TOF23/3600; % (hours) time of flight from SOI to periselenium
TOFhr=TOF01hr+TOF23hr; % (hours) total time of flight
altM(i,j)=RpM(i,j)-RM; % (SI) altitude at periselenium
vpM(i,j)=sqrt(2*(EM+muM/RpM(i,j))); % (SI) velocity at periselenium
if altM(i,j)>=999.9 && altM(i,j)<=1000.1 % altitude condition
    A(i,j)=1; % 1 if the condition is met, 0 if not
end
if vpM(i,j)<2.1 % velocity condition
    B(i,j)=1;
end
if A(i,j)==1 && B(i,j)==1
    iworks=i; % find interation that results in a successful lunar encounter
    jworks=j;
end
    end
end
HT=R0*v0(iworks); % (canonical) same parameters as before, but to check the candidate trajectory
ET=v0(iworks)^2/2-muE/R0; % (canonical)
pT=HT^2/muE; % (canonical)
aT=-muE/(2*ET); % (canonical)
eT=sqrt(1-pT/aT); % (dimensionless)
R1=sqrt(RsoiM^2+D^2-2*RsoiM*D*cosd(lambda1(jworks))); % (canonical)
v1=sqrt(2*(ET+muE/R1)); % (canonical)
phi1=acosd(real(HT/R1/v1)); % (degrees)
gamma1=asind(RsoiM/R1*sind(lambda1(jworks))); % (degrees)
nu0=real(acosd((pT-R0)/(R0*eT))); % (degrees)
nu1=real(acosd((pT-R1)/(R1*eT))); % (degrees)
EA0=acos((eT+cosd(nu0))/(1+eT*cosd(nu0))); % (radians)
EA1=acos((eT+cosd(nu1))/(1+eT*cosd(nu1))); % (radians)
TOF01=sqrt(aT^3/muE)*(EA1-eT*sin(EA1)-EA0+eT*sin(EA0)); % (canonical)
TOF01hr=TOF01*TU/3600; % (hours)
beta01=wM*TOF01; % (degrees)
gamma0=nu1-nu0-gamma1-wM*TOF01; % (degrees)
R2=RsoiMSI; % (SI)
vM=wMSI*DSI; % (SI)
v1SI=v1*DU/TU; % (SI)
v2=sqrt(v1SI^2+vM^2-2*v1SI*vM*cosd(real(phi1-gamma1))); % (SI)
eps2=asind(real(vM/v2*cosd(lambda1(jworks))-v1SI/v2*cosd(real(lambda1(jworks)+gamma1-phi1)))); % (degrees)
phi2=90-eps2; % (degrees)
EM=v2^2/2-muM/R2; % (SI)
HM=R2*v2*sind(real(eps2)); % (SI)
pM=HM^2/muM; % (SI)
eM=sqrt(1+2*EM*HM^2/muM^2); % (dimensionless)
aM=pM/(1-eM^2); % (SI)
RpM=pM/(1+eM); % (SI)
nu2=real(acosd(real((pM-R2)/(R2*eM)))); % (degrees)
nu3=real(acosd(real((pM-RpM)/(RpM*eM)))); % (degrees)
EA2=acos((eM+cosd(nu2))/(1+eM*cosd(nu2))); % (radians)
EA3=acos((eM+cosd(nu3))/(1+eM*cosd(nu3))); % (radians)
TOF23=-sqrt(aM^3/muM)*(EA3-eM*sin(EA3)-EA2+eM*sin(EA2)); % (SI)
TOF23hr=TOF23/3600; % (hours)
TOFhr=TOF01hr+TOF23hr; % (hours)
altMcheck=RpM-RM; % (SI)
vpMcheck=sqrt(2*(EM+muM/RpM)); % (SI)
if altMcheck==altM(iworks,jworks) % checks to see if candidate trajectory actually gets us there
    disp('Check')
    if vpMcheck==vpM(iworks,jworks)
        disp('Double Check')
        fprintf('\n %19s %.15f %s \n','v0 =',v0(iworks),'DU/TU') % read out parameters
        fprintf('%20s %.15f %s \n','phi0 =',phi0,'degrees')
        fprintf('%20s %.15f %s \n','lambda1 =',lambda1(jworks),'degrees')
        fprintf('%20s %.15f %s \n','Em0 =',ET,'DU^2/TU^2')
        fprintf('%20s %.15f %s \n','nu0 =',nu0,'degrees')
        fprintf('%20s %.15f %s \n','nu1 =',nu1,'degrees')
        fprintf('%20s %.15f %s \n','beta01 =',beta01,'degrees')
        fprintf('%20s %.15f %s \n','gamma0 =',gamma0,'degrees')
        fprintf('%20s %.15f %s \n','r1 =',R1,'DU')
        fprintf('%20s %.15f %s \n','v1 =',v1,'DU')
        fprintf('%20s %.15f %s \n','phi1 =',phi1,'degrees')
        fprintf('%20s %.15f %s \n','gamma1 =',gamma1,'degrees')
        fprintf('%20s %.15f %s \n','aCL =',aT,'DU')
        fprintf('%20s %.15f \n','eCL =',eT)
        fprintf('%20s %.15f %s \n','v2 =',v2,'km/s')
        fprintf('%20s %.15f %s \n','phi2 =',phi2,'degrees')
        fprintf('%20s %.15f %s \n','eps2 =',eps2,'degrees')
        fprintf('%20s %.15f %s \n','Em2 =',EM,'km^2/s^2')
        fprintf('%20s %.15f %s \n','a2 =',aM,'km')
        fprintf('%20s %.15f \n','e2 =',eM)
        fprintf('%20s %.15f %s \n','rp3 =',RpM,'km')
        fprintf('%20s %.15f %s \n','vp3 =',vpMcheck,'km/s')
        fprintf('%20s %.0f %s \n','TOF Earth to SOI =',TOF01hr,'hrs')
        fprintf('%20s %.0f %s \n','TOF SOI to Moon =',TOF23hr,'hrs')
        fprintf('%20s %.0f %s \n','Total TOF =',TOFhr,'hrs')
    else
        disp('No Check')
    end
else
    disp('No Check')
end