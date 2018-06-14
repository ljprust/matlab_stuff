flag=input('Canonical or km-s? (1 2) ');
if flag==1
    mu=1;
else
    mu=398600;
end

rmag=sqrt(r(1)^2+r(2)^2+r(3)^2);
vmag=sqrt(v(1)^2+v(2)^2+v(3)^2);

I=[1 0 0];
J=[0 1 0];
K=[0 0 1];

Hm=cross(r,v);
Hmmag=sqrt(Hm(1)^2+Hm(2)^2+Hm(3)^2);

Em=vmag^2/2-mu/rmag;

a=-mu/2/Em;

p=Hmmag^2/mu;

e=1/mu*((vmag^2-mu/rmag)*r-dot(r,v)*v);
emag=sqrt(1-p/a);
emag2=sqrt(e(1)^2+e(2)^2+e(3)^2);

n=cross(K,Hm);
nmag=sqrt(n(1)^2+n(2)^2+n(3)^2);

nu=acosd(dot(e,r)/emag/rmag);
if dot(r,v)<0
    nu=-nu;
end

inc=acosd(dot(Hm,K)/Hmmag);

OMEGA=acosd(dot(n,I)/nmag);
if n(2)<0
    OMEGA=-OMEGA;
end

omega=acosd(dot(n,e)/nmag/emag);
if e(3)<0
    omega=-omega;
end