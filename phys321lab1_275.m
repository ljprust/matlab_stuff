% Logan Prust
% Phys 321L Lab 1
% 275 Counts per 6 Seconds

clear,clc;

data=dlmread('phys321lab1_275.txt');
min=min(data);
max=max(data);
bins=linspace(min,max,10);
figure(1);
hist(data,bins);
title('Data (275 Counts per 6 Seconds)');
xlabel('Radiation Counts');
ylabel('Number of Occurrences');

ntrials=100;
j=1;
i=1;
mu=274.4;
sigma=17.35;
for k=min:max
    P(k)=mu^k*exp(-mu)/factorial(k)*ntrials;
    G(k)=1/sqrt(2*pi*sigma^2)*exp(-(k-mu)^2/2/sigma^2)*ntrials;
end

for k=min:max
    for m=1:P(k)
        fillP(j)=k;
        j=j+1;
    end
end

figure(2);
hist(fillP,bins);
title('Poisson Distribution (275 Counts per 6 Seconds)');
xlabel('Radiation Counts');
ylabel('Number of Occurrences');

for k=min:max
    for q=1:G(k)
        fillG(i)=k;
        i=i+1;
    end
end

figure(3);
hist(fillG,bins);
title('Gaussian Distribution (275 Counts per 6 Seconds)');
xlabel('Radiation Counts');
ylabel('Number of Occurrences');