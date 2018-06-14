% Logan Prust Aer E 261 Project 3 Main Program
clear,clc; % clears the workspace and screen
input=dlmread('alphaclinput.txt',' ',1,0); % inputs data points from text file
a=input(1,:); % extracts data point a
b=input(2,:); % extracts data point b
c=input(3,:); % extracts data point c
d=input(4,:); % extracts data point d
e=input(5,:); % extracts data point e
f=input(6,:); % extracts data point f
g=input(7,:); % extracts data point g
alpha=-g(1):.5:g(1); % sets up angle of attack matrix
for i=a(1):.5:b(1) % for-loop for first linear segment
    Cl(2*i+1)=alphacllinear(i,a(1),a(2),b(1),b(2)); % calculates Cl values in first linear segment
end % terminates for-loop
for i=b(1):.5:d(1) % for-loop for first quadratic segment
    Cl(2*i+1)=alphaclquadratic(i,b(1),b(2),c(1),c(2),d(1),d(2)); % calculates Cl values in first quadratic segment
end % terminates for-loop
for i=d(1):.5:f(1) % for-loop for second quadratic segment
    Cl(2*i+1)=alphaclquadratic(i,d(1),d(2),e(1),e(2),f(1),f(2)); % calculates Cl values in second quadratic segment
end % terminates for-loop
for i=f(1):.5:g(1) % for-loop for second linear segment
    Cl(2*i+1)=alphacllinear(i,f(1),f(2),g(1),g(2)); % calculates Cl values in second linear segment
end % terminates for-loop
negativeCl=fliplr(Cl(2:length(Cl))); % creates Cl matrix corresponding to negative angle of attack values
allCl=[negativeCl,Cl]; % synthesizes Cl matrices
plot(alpha,allCl); % plots Cl with respect to angle of attack
title('Angle of Attack vs. Coefficient of Lift'); % puts a title on the graph
xlabel('Angle of Attack (degrees)'); % puts a label on the x-axis
ylabel('Coefficient of Lift'); % puts a label on the y-axis