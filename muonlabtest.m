clear;

x=[1;2;3;4];
y=[1;3;4;4];
W=diag([1,0,1,0]);

n=length(x);

for i=1:n
    X(i,1)=x(i);
    X(i,2)=1;
end

betaW=inv(transpose(X)*W*X)*transpose(X)*W*y

beta=inv(transpose(X)*X)*transpose(X)*y;