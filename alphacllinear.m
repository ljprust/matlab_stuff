function[y]=alphacllinear(x,x1,y1,x2,y2)
y=(y2-y1)/(x2-x1)*x-(y2-y1)/(x2-x1)*x1+y1;
end