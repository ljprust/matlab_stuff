clear,clc;

A=[-2 3 20 0;0 3 -1 9;0 0 -6 1;0 0 0 7];
b=[7;8;9;10];
n=4;
C=[A b];

for i=1:n
    C(i,:)=C(i,:)/C(i,i);
end

for i=n:-1:2
    for j=1:i-1
        C(j,:)=-C(i,:)*C(j,i)+C(j,:);
    end
end