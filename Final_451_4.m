% Logan Prust - Aer E 451 - Problem Set Final - Problem 4

clear,clc;
mu=3.986*10^5;
tau1=-900;
tau3=900;
r2mag=13570;
v2mag=9.95;
r2=[1098,9974,9136];
v2=[-9.1,3.87,-1.01];
f1(1)=.9354;
f3(1)=.9354;
g1(1)=-880.6;
g3(1)=880.6;
for n=2:5
zeta=2/r2mag-v2mag^2/mu;
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
    if imag(rootcandidates1(i))==0
        x1=rootcandidates1(i);
    end
end
for i=1:13
    if imag(rootcandidates3(i))==0
        x3=rootcandidates3(i);
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
f1(n)=1-x1^2/r2mag*C1;
f3(n)=1-x3^2/r2mag*C1;
g1(n)=tau1-x1^3/sqrt(mu)*S1;
g3(n)=tau3-x3^3/sqrt(mu)*S3;
f1avg(n)=(f1(n-1)+f1(n))/2;
f3avg(n)=(f3(n-1)+f3(n))/2;
g1avg(n)=(g1(n-1)+g1(n))/2;
g3avg(n)=(g3(n-1)+g3(n))/2;
c1(n)=g3avg(n)/(f1avg(n)*g3avg(n)-f3avg(n)*g1avg(n));
c2(n)=g1avg(n)/(f1avg(n)*g3avg(n)-f3avg(n)*g1avg(n));
rho1(n)=1/D0*(-D(1,1)+D(2,1)/c1(n)-c3(n)/c1(n)*D(3,1));
rho2(n)=1/D0*(-c1(n)*D(1,2)+D(2,2)-c3(n)*D(3,2));
rho3(n)=1/D0*(-c1(n)/c3(n)*D(1,3)+D(2,3)/c3(n)-D(3,3));
r1=R(1,:)+rho1(n)*rhohat(1,:);
r2=R(2,:)+rho2(n)*rhohat(2,:);
r3=R(3,:)+rho3(n)*rhohat(3,:);
v2=(f1(n)*r3-f3(n)*r1)/(f1(n)*g3(n)-f3(n)*g3(n));
end