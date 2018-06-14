epsilon=10^(-6); % sets convergence criterion
x_initial=[0.001,2]; % sets initial values to start search
gamma=1.4;
% Astar=0.5009;
A=input('A/A* = ');
for j=1:2
    n=1; % initializes the iteration counter
    error=1; % sets some arbitrary initial value for the error
    x=x_initial(j); % sets x equal to one of its initial values
while error>=epsilon
    n=n+1; % increases the iteration counter
    f=1/x*(2/(gamma+1)*(1+(gamma-1)/2*x^2))^((gamma+1)/2/(gamma-1))-A; % evaluates f(x)
    fprime=-1/(x^2)*(2/(gamma+1)*(1+(gamma-1)/2*x^2))^((gamma+1)/2/(gamma-1))+(gamma+1)/2/(gamma-1)*1/x*(2/(gamma+1)*(1+(gamma-1)/2*x^2))^((gamma+1)/2/(gamma-1)-1)*2*(gamma-1)/(gamma+1)*x;
    % ^ evaluates f'(x)
    x2=x-f/fprime; % evaluates x(n+1)
    error=abs(x2-x); % evaluates the absolute error
    x=x2; % makes x(n+1) the new x(n)
end
fprintf('%20s %.15f \n','M =',x); % outputs the final value of x (also known as M)
xvector(j)=x; % puts the two M values found into vector form
end