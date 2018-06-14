% Logan Prust - Aer E 451 - Problem Set 2 - Problem 3

clear,clc;
x=zeros(1,4);
z=zeros(1,4);
C=zeros(1,4);
S=zeros(1,4);
F=zeros(1,4);
Fprime=zeros(1,4);
Fprimeprime=zeros(1,4);
delta=zeros(1,4);
signFprime=zeros(1,4);
mu=3.986e5; % km^3/s^2
fprintf('%20s %.15f %s \n','mu =',mu,'km^3/s^2');
t=60*60; % seconds
fprintf('%20s %.15f %s \n','t =',t,'seconds');
r0=[7000.0,-12124,0];
v0=[2.6679,4.6210,0];
r0mag=sqrt((r0(1))^2+(r0(2))^2+(r0(3))^2);
fprintf('%20s %.15f %s \n','r0 =',r0mag,'km');
v0mag=sqrt((v0(1))^2+(v0(2))^2+(v0(3))^2);
fprintf('%20s %.15f %s \n','v0 =',v0mag,'km/s');
a=mu/(2*mu/r0mag-v0mag^2);
fprintf('%20s %.15f %s \n','a =',a,'km');
x(1)=sqrt(mu)*t/abs(a);
for i=1:4
    z(i)=x(i)^2/a;
    C(i)=(1-cos(sqrt(z(i))))/z(i);
    S(i)=(sqrt(z(i))-sin(sqrt(z(i))))/(z(i))^(3/2);
    F(i)=(1-r0mag/a)*S(i)*(x(i))^3+dot(r0,v0)/sqrt(mu)*C(i)*(x(i))^2+r0mag*x(i)-sqrt(mu)*t;
    Fprime(i)=C(i)*(x(i))^2+dot(r0,v0)/sqrt(mu)*(1-S(i)*z(i))*x(i)+r0mag*(1-C(i)*z(i));
    Fprimeprime(i)=(1-r0mag/a)*(1-S(i)*z(i))*x(i)+dot(r0,v0)/sqrt(mu)*(1-C(i)*z(i));
    if Fprime<0
        signFprime(i)=-1;
    else
        signFprime(i)=1;
    end
    delta(i)=2*sqrt(4*(Fprime(i))^2-5*Fprime(i)*Fprimeprime(i));
    x(i+1)=x(i)-5*F(i)/(Fprime(i)+signFprime(i)*delta(i));
end
fprintf('%20s %.15f %s \n','x0 =',x(1),'km^(1/2)');
fprintf('%20s %.15f %s \n','x1 =',x(2),'km^(1/2)');
fprintf('%20s %.15f %s \n','x2 =',x(3),'km^(1/2)');
fprintf('%20s %.15f %s \n','x3 =',x(4),'km^(1/2)');
f=1-(x(4))^2*C(4)/r0mag;
fprintf('%20s %.15f %s \n','f =',f,'');
g=t-(x(4))^3*S(4)/sqrt(mu);
fprintf('%20s %.15f %s \n','g =',g,'seconds');
r=f*r0+g*v0;
fprintf('%20s %.15f %s \n','r_I =',r(1),'km');
fprintf('%20s %.15f %s \n','r_J =',r(2),'km');
fprintf('%20s %.15f %s \n','r_K =',r(3),'km');
rmag=sqrt((r(1))^2+(r(2))^2+(r(3))^2);
fprintf('%20s %.15f %s \n','r =',rmag,'km');
fdot=sqrt(mu)/rmag/r0mag*(x(4)*(z(4)*S(4)-1));
fprintf('%20s %.15f %s \n','f dot =',fdot,'seconds^-1');
gdot=1-(x(4))^2*C(4)/rmag;
fprintf('%20s %.15f %s \n','g dot =',gdot,'');
v=fdot*r0+gdot*v0;
fprintf('%20s %.15f %s \n','v_I =',v(1),'km/s');
fprintf('%20s %.15f %s \n','v_J =',v(2),'km/s');
fprintf('%20s %.15f %s \n','v_K =',v(3),'km/s');