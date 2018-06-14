% Logan Prust - Aer E 451 - Problem Set 4 - Problem 1

clear,clc;
digits(1000);
t=[0,60,120]; % sec
theta=[0,.250684,.501369]; % deg
RA=[0,65.9279,79.850]; % deg
declination=[51.5110,27.9911,14.6609]; % deg
phi=29; % deg
fprintf('%20s %.15f %s \n','phi (Reference Ellipsoid) =',phi,'deg');
tau1=t(2)-t(1); % sec
fprintf('%20s %.15f %s \n','tau1 =',tau1,'sec');
tau3=t(3)-t(2); % sec
fprintf('%20s %.15f %s \n','tau3 =',tau3,'sec');
tau=t(3)-t(1); % sec
fprintf('%20s %.15f %s \n','tau =',tau,'sec');
mu=3.986*10^5; % km^3/sec^2
fprintf('%20s %.15f %s \n','mu =',mu,'km^3/sec^2');
aRE=6378.1; % km
fprintf('%20s %.15f %s \n','a (Reference Ellipsoid) =',aRE,'km');
eRE=.08182;
fprintf('%20s %.15f %s \n','e (Reference Ellipsoid) =',eRE,'');
HRE=0; % km
fprintf('%20s %.15f %s \n','H (station altitude) =',HRE,'km');
xRE=(aRE/sqrt(1-eRE^2*sind(phi)^2)+HRE)*cosd(phi); % km
fprintf('%20s %.15f %s \n','x (Reference Ellipsoid =',xRE,'km');
zRE=(aRE*(1-eRE^2)/sqrt(1-eRE^2*sind(phi)^2)+HRE)*sind(phi); % km
fprintf('%20s %.15f %s \n','z (Reference Ellipsoid =',zRE,'km');
for i=1:3
    rhohat(i,1)=cosd(RA(i))*cosd(declination(i));
    rhohat(i,2)=sind(RA(i))*cosd(declination(i));
    rhohat(i,3)=sind(declination(i));
    R(i,1)=xRE*cosd(theta(i));
    R(i,2)=xRE*sind(theta(i));
    R(i,3)=zRE;
end
D0=dot(rhohat(1,:),cross(rhohat(2,:),rhohat(3,:)));
fprintf('%20s %.15f %s \n','D0 =',D0,'');
for i=1:3
    D(i,1)=dot(R(i,:),cross(rhohat(2,:),rhohat(3,:)));
    D(i,2)=dot(R(i,:),cross(rhohat(1,:),rhohat(3,:)));
    D(i,3)=dot(R(i,:),cross(rhohat(1,:),rhohat(2,:)));
end
A=1/D0*(-D(1,2)*tau3/tau+D(2,2)+D(3,2)*tau1/tau);
fprintf('%20s %.15f %s \n','A =',A,'km');
B=1/6/D0*(D(1,2)*(tau3^2-tau^2)*tau3/tau+D(3,2)*(tau^2-tau1^2)*tau1/tau);
fprintf('%20s %.15f %s \n','B =',B,'km*sec^2');
E=dot(R(2,:),rhohat(2,:));
fprintf('%20s %.15f %s \n','E =',E,'km');
R2=sqrt(R(2,1)^2+R(2,2)^2+R(2,3)^2);
acoeff=-(A^2+2*A*E+R2^2);
fprintf('%20s %.15f %s \n','a (coefficient) =',acoeff,'km^2');
bcoeff=-2*mu*B*(A+E);
fprintf('%20s %.15f %s \n','b (coefficient) =',bcoeff,'km^5');
ccoeff=-mu^2*B^2;
fprintf('%20s %.15f %s \n','c (coefficient) =',ccoeff,'km^8');
mayber2=roots([1,0,acoeff,0,0,bcoeff,0,0,ccoeff]);
for i=1:8
    if imag(mayber2(i))==0
        if mayber2(i)>=0
            r2mag=mayber2(i);
        end
    end
end
fprintf('%20s %.15f %s \n','r2 =',r2mag,'km');
f1=1-mu/2/r2mag^3*tau1^2;
fprintf('%20s %.15f %s \n','f1 =',f1,'');
f3=1-mu/2/r2mag^3*tau3^2;
fprintf('%20s %.15f %s \n','f3 =',f3,'');
g1=tau1-mu/6/r2mag^3*tau1^3;
fprintf('%20s %.15f %s \n','g1 =',g1,'sec');
g3=tau3-mu/6/r2mag^3*tau3^3;
fprintf('%20s %.15f %s \n','g3 =',g3,'sec');
Ldenominator=tau-mu/6/r2mag^3*tau^3;
fprintf('%20s %.15f %s \n','Lagrange denominator =',Ldenominator,'sec');
c1=g3/Ldenominator;
fprintf('%20s %.15f %s \n','c1 =',c1,'');
c3=-g1/Ldenominator;
fprintf('%20s %.15f %s \n','c3 =',c3,'');
rho1=1/D0*(-D(1,1)+D(2,1)/c1-c3/c1*D(3,1));
fprintf('%20s %.15f %s \n','rho1 =',rho1,'km');
rho2=1/D0*(-c1*D(1,2)+D(2,2)-c3*D(3,2));
fprintf('%20s %.15f %s \n','rho2 =',rho2,'km');
rho3=1/D0*(-c1/c3*D(1,3)+D(2,3)/c3-D(3,3));
fprintf('%20s %.15f %s \n','rho3 =',rho3,'km');
r1=R(1,:)+rho1*rhohat(1,:);
r1mag=sqrt(r1(1)^2+r1(2)^2+r1(3)^2);
fprintf('%20s %.15f %s \n','r1 =',r1mag,'km');
r3=R(3,:)+rho3*rhohat(3,:);
r3mag=sqrt(r3(1)^2+r3(2)^2+r3(3)^2);
fprintf('%20s %.15f %s \n','r3 =',r3mag,'km');
v2=(f1*r3-f3*r1)/Ldenominator;
v2mag=sqrt(v2(1)^2+v2(2)^2+v2(3)^2);
fprintf('%20s %.15f %s \n','v2 =',v2mag,'km/sec');
fprintf('%20s %.15f %s \n','v2_I =',v2(1),'km/sec');
fprintf('%20s %.15f %s \n','v2_J =',v2(2),'km/sec');
fprintf('%20s %.15f %s \n','v2_K =',v2(3),'km/sec');
r2=1/f1*r1-g1/f1*v2;
fprintf('%20s %.15f %s \n','r2_I =',r2(1),'km');
fprintf('%20s %.15f %s \n','r2_J =',r2(2),'km');
fprintf('%20s %.15f %s \n','r2_K =',r2(3),'km');