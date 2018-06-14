e=0.2;
p=1-e^2;
nu=linspace(0,2*pi,100);

for i=1:100
    r(i)=p/(1+e*cos(nu(i)));
end

polar(nu,r);