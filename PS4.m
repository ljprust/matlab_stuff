% Aer E 351 PS 4 #1
clear,clc;
p=.125;
e=.8839;
nu=1:360;
r=zeros(1,360);
for i=1:360
    r(i)=p/(1+e*cosd(nu(i)));
end
polar(nu*3.14/180,r);