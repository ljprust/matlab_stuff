x=[.12,.18,.25,.32,.38];
y=[.51,.73,.85,.79,1.14];
sigmay=[.07,.11,.15,.18,.12];
a1=(5*dot(x,y)-x(1)*x(2)*x(3)*x(4)*x(5))/(5*x(1)^2*x(2)^2*x(3)^2*x(4)^2*x(5)^2-(x(1)*x(2)*x(3)*x(4)*x(5))^2)
b1=((x(1)*x(2)*x(3)*x(4)*x(5))^2*y(1)*y(2)*y(3)*y(4)*y(5)-x(1)*x(2)*x(3)*x(4)*x(5)*dot(x,y))/(5*x(1)^2*x(2)^2*x(3)^2*x(4)^2*x(5)^2-(x(1)*x(2)*x(3)*x(4)*x(5))^2)
s=0;
for i=1:5
    s=s+(y(i)-a1*x(i)-b1)^2;
end
S=sqrt(s/3);
sigmaa1=S*sqrt(5/(5*x(1)^2*x(2)^2*x(3)^2*x(4)^2*x(5)^2-(x(1)*x(2)*x(3)*x(4)*x(5))^2))
sigmab1=S*sqrt((x(1)^2*x(2)^2*x(3)^2*x(4)^2*x(5)^2)/(5*x(1)^2*x(2)^2*x(3)^2*x(4)^2*x(5)^2-(x(1)*x(2)*x(3)*x(4)*x(5))^2))
for j=1:5
    blah1=x(j)/sigmay(j)^2;
    blah2=y(j)/sigmay(j)^2;
    blah3=x(j)*y(j)/sigmay(j)^2;
    blah4=1/sigmay(j)^2;
    blah5=x(j)^2/sigmay(j)^2;
end
a2=(blah1*blah2-blah3*blah4)/(blah1^2-blah5*blah4)
b2=(blah1*blah3-blah2*blah5)/(blah1^2-blah5*blah4)
sigmaa2=sqrt(blah4/(blah5*blah4-blah1^2))
sigmab2=sqrt(blah5/(blah5*blah4-blah1^2))
disp(blah1)
disp(blah2)
disp(blah3)
disp(blah4)
disp(blah5)