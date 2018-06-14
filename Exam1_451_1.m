% Logan Prust - Aer E 451 - Exam 1 - Problem 1

clear,clc;
mu=3.985*10^5; % km^3/sec^2
TOF=9141.513692; % sec
r1=[0,30000,0]; % km
r1mag=30000; % km
r2=[-22500,0,0]; % km
r2mag=22500; % km
deltanu=90; % deg
A=sqrt(r1mag*r2mag)*sind(deltanu)/sqrt(1-cosd(deltanu));
fprintf('%20s %.15f %s \n','A =',A,'km');
n=1000000;
z=linspace(0,pi^2,n);
for i=1:n
    C(i)=(1-cos(sqrt(z(i))))/z(i);
    S(i)=(sqrt(z(i))-sin(sqrt(z(i))))*(z(i))^(-3/2);
    y(i)=r1mag+r2mag-A*(1-z(i)*S(i))/sqrt(C(i));
    x(i)=sqrt(y(i)/C(i));
    t(i)=(x(i))^3*S(i)/sqrt(mu)+A*sqrt(y(i))/sqrt(mu);
    if t(i)>TOF-.01 && t(i)<TOF+.01
        ifinal=i;
    end
end
semilogy(z,t);
fprintf('%20s %.15f %s \n','z final =',z(ifinal),'');
f=1-y(ifinal)/r1mag;
fprintf('%20s %.15f %s \n','f =',f,'');
g=A*sqrt(y(ifinal)/mu);
fprintf('%20s %.15f %s \n','g =',g,'sec');
gdot=1-y(ifinal)/r2mag;
fprintf('%20s %.15f %s \n','g dot =',gdot,'');
v1=(r2-f*r1)/g;
fprintf('%20s %.15f %s \n','v1P =',v1(1),'km/sec');
fprintf('%20s %.15f %s \n','v1Q =',v1(2),'km/sec');
fprintf('%20s %.15f %s \n','v1W =',v1(3),'km/sec');
v1mag=sqrt(v1(1)^2+v1(2)^2+v1(3)^2);
fprintf('%20s %.15f %s \n','v1 =',v1mag,'km/sec');
v2=(gdot*r2-r1)/g;
fprintf('%20s %.15f %s \n','v2P =',v2(1),'km/sec');
fprintf('%20s %.15f %s \n','v2Q =',v2(2),'km/sec');
fprintf('%20s %.15f %s \n','v2W =',v2(3),'km/sec');
v2mag=sqrt(v2(1)^2+v2(2)^2+v2(3)^2);
fprintf('%20s %.15f %s \n','v2 =',v2mag,'km/sec');