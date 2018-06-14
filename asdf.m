nu = linspace(0,2*pi,1000);
e = 0.2;
p = 0.3;
for i = 1:1000
    r(i) = p/(1+e*cos(nu(i)));
end

polar(nu,r);