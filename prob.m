% Probability
clear,clc;
results=zeros(1,12);
iterations=1000000;
n=zeros(1,iterations);
for i=1:iterations
    n(i)=randi(6)+randi(6);
    results(n(i))=results(n(i))+1;
end
output=results'/iterations*100;
disp(output);
bounds=1:12;
hist(n,bounds);