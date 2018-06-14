% Logan Prust
% Phys 321L Lab 1
% Penny Tossing

clear,clc;

pennies=dlmread('Penny Counting.txt');
bins=0:7;
figure(1);
hist(pennies,bins);
title('Penny Data');
xlabel('Number of Heads');
ylabel('Number of Occurrences');

ntrials=200;
j=1;
n=7;
p=0.5;
for k=0:7
    B(k+1)=factorial(n)/factorial(k)/factorial(n-k)*p^n*ntrials;
end

for k=0:7
    for m=1:B(k+1)
        fill(j)=k;
        j=j+1;
    end
end

figure(2);
hist(fill,bins);
title('Binomial Distribution for Pennies');
xlabel('Number of Heads');
ylabel('Number of Occurrences');