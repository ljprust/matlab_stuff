% Logan Prust - Aer E 451 - Exam 1 - Problem 2

clear,clc;
inclination=28.5; % deg
mu=3.986*10^5; % km^3/sec^2
m=100000;
alpha1=linspace(0,inclination,m); % deg
r0=6700; % km
rGSO=42164; % km
a=(r0+rGSO)/2; % km
fprintf('%20s %.15f %s \n','a =',a,'km');
vT0=sqrt(2*mu/r0-mu/a); % km/sec
fprintf('%20s %.15f %s \n','vT0 =',vT0,'km/sec');
vTGSO=sqrt(2*mu/rGSO-mu/a); % km/sec
fprintf('%20s %.15f %s \n','vTGSO =',vTGSO,'km/sec');
v0=sqrt(mu/r0); % km/sec
fprintf('%20s %.15f %s \n','vcirc0 =',v0,'km/sec');
vGSO=sqrt(mu/rGSO); % km/sec
fprintf('%20s %.15f %s \n','vcircGSO =',vGSO,'km/sec');
ToF=pi*sqrt(a^3/mu)/3600; % hrs
fprintf('%20s %.15f %s \n','ToF =',ToF,'hours');
for i=1:m
    DV1(i)=sqrt(v0^2+vT0^2-2*v0*vT0*cosd(alpha1(i)));
    DV2(i)=sqrt(vTGSO^2+vGSO^2-2*vTGSO*vGSO*cosd(inclination-alpha1(i)));
    DV(i)=DV1(i)+DV2(i);
    if i==1
        imin=1;
    else
        if DV(i)<DV(i-1)
            imin=i;
        end
    end
end
C=DV(imin)+2.0*ToF;
fprintf('%20s %.15f %s \n','alpha 1 min =',alpha1(imin),'degrees');
fprintf('%20s %.15f %s \n','delta v min =',DV(imin),'km/sec');
fprintf('%20s %.15f %s \n\n','C min =',C,'km');

clear;
inclination=28.5; % deg
mu=3.986*10^5; % km^3/sec^2
m=100;
alpha1=linspace(0,inclination,m); % deg
alpha2=inclination-alpha1; % den
r0=6700; % km
rGSO=42164; % km
v0=sqrt(mu/r0); % km/sec
vGSO=[0,sqrt(mu/rGSO),0]; % km/sec
aHohmann=(r0+rGSO)/2; % km
vT0Hohmann=sqrt(2*mu/r0-mu/aHohmann);
n=100;
vT0=linspace(vT0Hohmann,10.9,n); % km/sec
Cmin=10000000;
for i=1:m
    for j=1:n
        a(j)=mu/(2*mu/r0-vT0(j)^2);
        p(j)=(r0^2*vT0(j)^2)/mu;
        e(j)=sqrt(1-p(j)/a(j));
        vTGSO(j)=sqrt(2*mu/rGSO-mu/a(j));
        phi(j)=acosd(r0*vT0(j)/rGSO/vTGSO(j));
        vx=vTGSO(j)*cosd(alpha2(i))*sind(phi(j));
        vy=vTGSO(j)*cosd(alpha2(i))*cosd(phi(j));
        vz=-vTGSO(j)*sind(alpha2(i));
        v=[vx,vy,vz];
        dv=vGSO-v;
        DV1(i,j)=sqrt(v0^2+vT0(j)^2-2*v0*vT0(j)*cosd(alpha1(i)));
        DV2(i,j)=sqrt(dv(1)^2+dv(2)^2+dv(3)^2);
        DV(i,j)=DV1(i,j)+DV2(i,j);
        nu(j)=acosd(p(j)/rGSO/e(j)-1/e(j));
        E(j)=2*atan(sqrt((1-e(j))/(1+e(j)))*tand(nu(j)/2));
        if E(j)<0
            E(j)=E(j)+pi/2;
        end
        ToF(j)=sqrt(a(j)^3/mu)*(E(j)-e(j)*sin(E(j)))/3600; % hrs
        C(i,j)=DV(i,j)+2.0*ToF(j);
        if C(i,j)<Cmin
            Cmin=C(i,j);
            imin=i;
            jmin=j;
        end
    end
end
fprintf('%20s %.15f %s \n','delta nu min =',nu(jmin),'degrees');
fprintf('%20s %.15f %s \n','delta v min =',DV(imin,jmin),'km/sec');
fprintf('%20s %.15f %s \n','alpha 1 min =',alpha1(imin),'degrees');
fprintf('%20s %.15f %s \n','alpha 2 min =',alpha2(imin),'degrees');
fprintf('%20s %.15f %s \n','ToF min =',ToF(jmin),'hours');
fprintf('%20s %.15f %s \n','C min =',C(imin,jmin),'');