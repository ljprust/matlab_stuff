% Logan Prust - Phys 321L - Lab 3

clear,clc;

nu(1)=3E8/365E-9;
nu(2)=3E8/405E-9;
nu(3)=3E8/436E-9;
nu(4)=3E8/546E-9;
nu(5)=3E8/577E-9;
nu(6)=nu(1);
nu(7)=nu(2);
nu(8)=nu(3);
nu(9)=nu(4);
nu(10)=nu(5);
nu(11)=nu(1);
nu(12)=nu(2);
nu(13)=nu(3);
nu(14)=nu(4);
nu(15)=nu(5);

V=[1.632 1.269 1.081 .578 .459 1.723 1.321 1.118 .598 .483 1.751 1.342 1.137 .605 .492];
V=-V;

sigmaV=0.0005;

n=15;
sum1=0;
sumx=0;
sumy=0;
sumx2=0;
sumxy=0;
for i=1:n
    sum1=sum1+1;
    sumx=sumx+nu(i);
    sumy=sumy+V(i);
    sumx2=sumx2+nu(i)^2;
    sumxy=sumxy+nu(i)*V(i);
end

a=(sumy*sumx2-sumx*sumxy)/(sum1*sumx2-sumx*sumx);
b=(sum1*sumxy-sumx*sumy)/(sum1*sumx2-sumx*sumx);
sigmaa2=sigmaV^2*sumx2/(sum1*sumx2-sumx*sumx);
sigmab2=sigmaV^2*sum1/(sum1*sumx2-sumx*sumx);