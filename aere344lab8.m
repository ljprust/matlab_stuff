rho=1.22;
mu=1.8*10^(-5);
x=input('x = ');
p_total=input('p total = ');
p_inf=input('p inf = ');
u_inf=sqrt(2/rho*(p_total-p_inf));
Re=rho*u_inf*x/mu;
cf=0.0583/Re^(1/5);
if Re<=5*10^5
    delta=5*x/sqrt(Re);
else
    delta=0.37*x/Re^(1/5);
end
disp('delta=')
disp(delta)
disp('cf=')
disp(cf)
u_delta=u_inf*0.99;
p_delta=p_inf+0.5*rho*u_delta^2;
disp('p_delta = ');
disp(p_delta);