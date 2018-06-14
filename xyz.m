x=7;
y=9;
z=38;
n=1;
seed=randi(x+y+z);
if seed<=x
    n=randi(x);
    fprintf('%5s %f \n','X',n);
elseif seed<=x+y
    n=randi(y);
    fprintf('%5s %f \n','Y',n);
else
    n=randi(z);
    fprintf('%5s %f \n','Z',n);
end