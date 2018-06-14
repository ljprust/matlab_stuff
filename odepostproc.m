% Logan Prust - Aer E 361 - ODE Post-Processor

clear,clc;

% get step size

H=input('Step size: ');

% get info from FORTRAN code

odeinput=dlmread('odedata.txt');
odeoutput=dlmread('odeoutput.txt');

% find number of state variables

N=odeinput(1);

% find initial conditions

for i=1:N
    init(i)=odeinput(i+1);
end

% find number of steps

M=odeinput(2+N);

% create array of x values

x=0:H:H*(M-1);

% construct arrays of y values for each state variable

k=0;
for i=1:N
    y(:,i)=odeoutput((1+k*M):(1000+k*M))';
    k=k+1;
end

% plot each state variable

for i=1:N
    figure(i);
    plot(x,y(:,i));
    xlabel('x');
    ylabel('y');
end

% plot all state variables on one plot

figure(i+1);
hold on
for j=1:N
    plot(x,y(:,j));
end
xlabel('x');
ylabel('y');