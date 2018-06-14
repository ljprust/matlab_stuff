% Logan Prust - Phys 481 - Homework 1 - Problem 5
clear,clc;
z=1;
error(1)=1;
while error(z)>=0.01
z=z+1;
error(z)=(log(factorial(z))-z*log(z)+z)/log(factorial(z));
end
fprintf('%.5f\n',z);