% Logan Prust - Aer E 451 - Problem Set Final - Problem 3

clear,clc;
fprintf('%20s \n\n','Gauss Method:');
digits(1000);
t=[0,15*60,30*60]; % sec
theta=[35.93863185,39.69889901,43.45916616]; % deg
RA=[55.0931551,111.4115235,134.2826693]; % deg
declination=[36.5731946,25.9983390,12.0351097]; % deg
phi=40; % deg
fprintf('%20s %.15f %s \n','phi (Reference Ellipsoid) =',phi,'deg');
tau1=t(1)-t(2); % sec
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
HRE=2; % km
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
fprintf('%20s %.5f %s \n','rhohat1 I =',rhohat(1,1),'');
fprintf('%20s %.5f %s \n','rhohat1 J =',rhohat(1,2),'');
fprintf('%20s %.5f %s \n','rhohat1 K =',rhohat(1,3),'');
fprintf('%20s %.5f %s \n','rhohat2 I =',rhohat(2,1),'');
fprintf('%20s %.5f %s \n','rhohat2 J =',rhohat(2,2),'');
fprintf('%20s %.5f %s \n','rhohat2 K =',rhohat(2,3),'');
fprintf('%20s %.5f %s \n','rhohat3 I =',rhohat(3,1),'');
fprintf('%20s %.5f %s \n','rhohat3 J =',rhohat(3,2),'');
fprintf('%20s %.5f %s \n','rhohat3 K =',rhohat(3,3),'');
fprintf('%20s %.5f %s \n','R1 I =',R(1,1),'km');
fprintf('%20s %.5f %s \n','R1 J =',R(1,2),'km');
fprintf('%20s %.5f %s \n','R1 K =',R(1,3),'km');
fprintf('%20s %.5f %s \n','R2 I =',R(2,1),'km');
fprintf('%20s %.5f %s \n','R2 J =',R(2,2),'km');
fprintf('%20s %.5f %s \n','R2 K =',R(2,3),'km');
fprintf('%20s %.5f %s \n','R3 I =',R(3,1),'km');
fprintf('%20s %.5f %s \n','R3 J =',R(3,2),'km');
fprintf('%20s %.5f %s \n','R3 K =',R(3,3),'km');
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
R1=sqrt(R(1,1)^2+R(1,2)^2+R(1,3)^2);
R2=sqrt(R(2,1)^2+R(2,2)^2+R(2,3)^2);
R3=sqrt(R(3,1)^2+R(3,2)^2+R(3,3)^2);
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
P=6*(D(3,1)*tau1/tau3+D(2,1)*tau/tau3)*r2mag^3+mu*D(3,1)*(tau^2-tau1^2)*tau1/tau3;
fprintf('%20s %.15f %s \n','P =',P,'km^4');
Q=6*(D(1,3)*tau3/tau1-D(2,3)*tau/tau1)*r2mag^3+mu*D(1,3)*(tau^2-tau3^2)*tau3/tau1;
fprintf('%20s %.15f %s \n','Q =',Q,'km^4');
scriptR1=6*r2mag^3+mu*(tau^2-tau3^2);
scriptR3=6*r2mag^3+mu*(tau^2-tau1^2);
rho1=1/D0*(P/scriptR1-D(1,1));
fprintf('%20s %.15f %s \n','rho1 =',rho1,'km');
rho2=A+mu/r2mag^3*B;
fprintf('%20s %.15f %s \n','rho2 =',rho2,'km');
rho3=1/D0*(Q/scriptR3-D(3,3));
fprintf('%20s %.15f %s \n','rho3 =',rho3,'km');
f1(1)=1-mu/2/r2mag^3*tau1^2;
fprintf('%20s %.15f %s \n','f1 =',f1,'');
f3(1)=1-mu/2/r2mag^3*tau3^2;
fprintf('%20s %.15f %s \n','f3 =',f3,'');
g1(1)=tau1-mu/6/r2mag^3*tau1^3;
fprintf('%20s %.15f %s \n','g1 =',g1,'sec');
g3(1)=tau3-mu/6/r2mag^3*tau3^3;
fprintf('%20s %.15f %s \n','g3 =',g3,'sec');
r1=R(1,:)+rho1*rhohat(1,:);
r1mag=sqrt(r1(1)^2+r1(2)^2+r1(3)^2);
fprintf('%20s %.15f %s \n','r1 =',r1mag,'km');
r2=R(2,:)+rho2*rhohat(2,:);
fprintf('%20s %.15f %s \n','r2 =',r2mag,'km');
r3=R(3,:)+rho3*rhohat(3,:);
r3mag=sqrt(r3(1)^2+r3(2)^2+r3(3)^2);
fprintf('%20s %.15f %s \n','r3 =',r3mag,'km');
v2=(f1*r3-f3*r1)/(f1*g3-f3*g1);
v2mag=sqrt(v2(1)^2+v2(2)^2+v2(3)^2);
fprintf('%20s %.15f %s \n','v2 =',v2mag,'km/sec');
fprintf('%20s %.15f %s \n','v2_I =',v2(1),'km/sec');
fprintf('%20s %.15f %s \n','v2_J =',v2(2),'km/sec');
fprintf('%20s %.15f %s \n','v2_K =',v2(3),'km/sec');
fprintf('\n %20s \n\n','Orbital Elements:');
Hm=cross(r2,v2);
Hmmag=sqrt(Hm(1)^2+Hm(2)^2+Hm(3)^2);
fprintf('%20s %.15f %s \n','Hm =',Hmmag,'km^2/sec^2');
nvector=cross([0,0,1],Hm);
nmag=sqrt(nvector(1)^2+nvector(2)^2+nvector(3)^2);
a=mu/(2*mu/r2mag-v2mag^2);
fprintf('%20s %.15f %s \n','a =',a,'km');
inc=acosd(dot(Hm,[0,0,1])/Hmmag);
fprintf('%20s %.15f %s \n','inclination =',inc,'deg');
Omega=acosd(dot([1,0,0],nvector)/nmag);
fprintf('%20s %.15f %s \n','Big Omega =',Omega,'deg');
evector=1/mu*((v2mag^2-mu/r2mag)*r2-dot(r2,v2)*v2);
e=sqrt(evector(1)^2+evector(2)^2+evector(3)^2);
fprintf('%20s %.15f %s \n','e =',e,'');
omega=acosd(dot(nvector,evector)/nmag/e);
fprintf('%20s %.15f %s \n','little omega =',omega,'deg');
nu0=acosd(dot(r2,evector)/r2mag/e);
fprintf('%20s %.15f %s \n','true anomaly angle =',nu0,'deg');

fprintf('\n %20s \n','Gauss Method With Improvement:');
for n=2:21
fprintf('\n %20s %1.f %s \n\n','Iteration',n-1,'of Improved Method');
zeta=2/r2mag-v2mag^2/mu;
fprintf('%20s %.15f %s \n','zeta =',zeta,'1/km');
a=1/zeta;
fprintf('%20s %.15f %s \n','a =',a,'km');
rootvector1(1)=-zeta^5/6227020800*(1-zeta*r2mag);
rootvector1(2)=-zeta^5/479001600*dot(r2,v2)/sqrt(mu);
rootvector1(3)=zeta^4/39916800*(1-zeta*r2mag);
rootvector1(4)=zeta^4/3628800*dot(r2,v2)/sqrt(mu);
rootvector1(5)=-zeta^3/362880*(1-zeta*r2mag);
rootvector1(6)=-zeta^3/40320*dot(r2,v2)/sqrt(mu);
rootvector1(7)=zeta^2/5040*(1-zeta*r2mag);
rootvector1(8)=zeta^2/720*dot(r2,v2)/sqrt(mu);
rootvector1(9)=-zeta/120*(1-zeta*r2mag);
rootvector1(10)=-zeta/24*dot(r2,v2)/sqrt(mu);
rootvector1(11)=1/6*(1-zeta*r2mag);
rootvector1(12)=1/2*dot(r2,v2)/sqrt(mu);
rootvector1(13)=r2mag;
rootvector1(14)=-1*tau1*sqrt(mu);
rootvector3(1)=-zeta^5/6227020800*(1-zeta*r2mag);
rootvector3(2)=-zeta^5/479001600*dot(r2,v2)/sqrt(mu);
rootvector3(3)=zeta^4/39916800*(1-zeta*r2mag);
rootvector3(4)=zeta^4/3628800*dot(r2,v2)/sqrt(mu);
rootvector3(5)=-zeta^3/362880*(1-zeta*r2mag);
rootvector3(6)=-zeta^3/40320*dot(r2,v2)/sqrt(mu);
rootvector3(7)=zeta^2/5040*(1-zeta*r2mag);
rootvector3(8)=zeta^2/720*dot(r2,v2)/sqrt(mu);
rootvector3(9)=-zeta/120*(1-zeta*r2mag);
rootvector3(10)=-zeta/24*dot(r2,v2)/sqrt(mu);
rootvector3(11)=1/6*(1-zeta*r2mag);
rootvector3(12)=1/2*dot(r2,v2)/sqrt(mu);
rootvector3(13)=r2mag;
rootvector3(14)=-1*tau3*sqrt(mu);
rootcandidates1=roots(rootvector1);
rootcandidates3=roots(rootvector3);
for i=1:13
    if imag(rootcandidates1(i))==0 && real(rootcandidates1(i))<0 && real(rootcandidates1(i))>-100
        x1=rootcandidates1(i);
        fprintf('%20s %.15f %s \n','x1 =',x1,'km^(1/2)');
    end
end
for i=1:13
    if imag(rootcandidates3(i))==0 && real(rootcandidates3(i))>0 && real(rootcandidates3(i))<100
        x3=rootcandidates3(i);
        fprintf('%20s %.15f %s \n','x3 =',x3,'km^(1/2)');
    end
end
if zeta>0
    C1=(1-cos(sqrt(zeta*x1^2)))/(zeta*x1^2);
    C3=(1-cos(sqrt(zeta*x3^2)))/(zeta*x3^2);
    S1=(sqrt(zeta*x1^2)-sin(sqrt(zeta*x1^2)))/(zeta*x1^2)^(3/2);
    S3=(sqrt(zeta*x3^2)-sin(sqrt(zeta*x3^2)))/(zeta*x3^2)^(3/2);
elseif zeta<0
    C1=(cosh(sqrt(-1*zeta*x1^2))-1)/(-1*zeta*x1^2);
    C3=(cosh(sqrt(-1*zeta*x3^2))-1)/(-1*zeta*x3^2);
    S1=(sinh(sqrt(-1*zeta*x1^2))-sqrt(-1*zeta*x1^2))/(-1*zeta*x1^2)^(3/2);
    S3=(sinh(sqrt(-1*zeta*x3^2))-sqrt(-1*zeta*x3^2))/(-1*zeta*x3^2)^(3/2);
else
    C1=1/2;
    C3=1/2;
    S1=1/6;
    S3=1/6;
end
fprintf('%20s %.15f %s \n','C1 (x1) =',C1,'');
fprintf('%20s %.15f %s \n','C3 (x3) =',C3,'');
fprintf('%20s %.15f %s \n','S1 (x1) =',S1,'');
fprintf('%20s %.15f %s \n','S3 (x3) =',S3,'');
f1(n)=1-x1^2/r2mag*C1;
fprintf('%20s %.15f %s \n','f1 =',f1(n),'');
f3(n)=1-x3^2/r2mag*C3;
fprintf('%20s %.15f %s \n','f3 =',f3(n),'');
g1(n)=tau1-x1^3/sqrt(mu)*S1;
fprintf('%20s %.15f %s \n','g1 =',g1(n),'sec');
g3(n)=tau3-x3^3/sqrt(mu)*S3;
fprintf('%20s %.15f %s \n','g3 =',g3(n),'sec');
f1avg(n)=(f1(n-1)+f1(n))/2;
f3avg(n)=(f3(n-1)+f3(n))/2;
g1avg(n)=(g1(n-1)+g1(n))/2;
g3avg(n)=(g3(n-1)+g3(n))/2;
c1(n)=g3avg(n)/(f1avg(n)*g3avg(n)-f3avg(n)*g1avg(n));
c3(n)=-1*g1avg(n)/(f1avg(n)*g3avg(n)-f3avg(n)*g1avg(n));
rho1(n)=1/D0*(-D(1,1)+D(2,1)/c1(n)-c3(n)/c1(n)*D(3,1));
rho2(n)=1/D0*(-c1(n)*D(1,2)+D(2,2)-c3(n)*D(3,2));
rho3(n)=1/D0*(-c1(n)/c3(n)*D(1,3)+D(2,3)/c3(n)-D(3,3));
r1=R(1,:)+rho1(n)*rhohat(1,:);
r1mag=sqrt(r1(1)^2+r1(2)^2+r1(3)^2);
fprintf('%20s %.15f %s \n','r1 =',r1mag,'km');
r2=R(2,:)+rho2(n)*rhohat(2,:);
r2mag=sqrt(r2(1)^2+r2(2)^2+r2(3)^2);
fprintf('%20s %.15f %s \n','r2 =',r2mag,'km');
r3=R(3,:)+rho3(n)*rhohat(3,:);
r3mag=sqrt(r3(1)^2+r3(2)^2+r3(3)^2);
fprintf('%20s %.15f %s \n','r3 =',r3mag,'km');
v2=(f1(n)*r3-f3(n)*r1)/(f1(n)*g3(n)-f3(n)*g1(n));
v2mag=sqrt(v2(1)^2+v2(2)^2+v2(3)^2);
fprintf('%20s %.15f %s \n','v2 =',v2mag,'km/sec');
end
fprintf('\n %20s \n\n','Final Orbital Elements:');
Hm=cross(r2,v2);
Hmmag=sqrt(Hm(1)^2+Hm(2)^2+Hm(3)^2);
fprintf('%20s %.15f %s \n','Hm =',Hmmag,'km^2/sec^2');
nvector=cross([0,0,1],Hm);
nmag=sqrt(nvector(1)^2+nvector(2)^2+nvector(3)^2);
a=mu/(2*mu/r2mag-v2mag^2);
fprintf('%20s %.15f %s \n','a =',a,'km');
inc=acosd(dot(Hm,[0,0,1])/Hmmag);
fprintf('%20s %.15f %s \n','inclination =',inc,'deg');
Omega=acosd(dot([1,0,0],nvector)/nmag);
fprintf('%20s %.15f %s \n','Big Omega =',Omega,'deg');
evector=1/mu*((v2mag^2-mu/r2mag)*r2-dot(r2,v2)*v2);
e=sqrt(evector(1)^2+evector(2)^2+evector(3)^2);
fprintf('%20s %.15f %s \n','e =',e,'');
omega=acosd(dot(nvector,evector)/nmag/e);
fprintf('%20s %.15f %s \n','little omega =',omega,'deg');
nu0=acosd(dot(r2,evector)/r2mag/e);
fprintf('%20s %.15f %s \n','true anomaly angle =',nu0,'deg');