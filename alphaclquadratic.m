function[y]=alphaclquadratic(x,x1,y1,x2,y2,x3,y3)
A=1/(x2-x3)*((y2-y1)/(x2-x1)-(y3-y1)/(x3-x1));
B=(y2-y1)/(x2-x1)-A*(x2+x1);
C=y1-A*x1^2-B*x1;
y=A*x^2+B*x+C;
end