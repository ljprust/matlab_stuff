% Logan Prust - Math 407 - HW 6

A=[4 8 2 1 -9 -1;1 -8 -2 -3 10 0;4 -15 -1 -6 15 -1;5 -10 -2 -2 10 -3;3 -15 -3 -5 17 -1;3 8 2 1 -9 0];

B=[1 1 1 1 1 1;0 1 2 3 4 5;1 2 4 8 16 32;0 1 4 12 32 80;0 0 2 12 48 160;0 0 0 6 48 240];

C=[exp(1);exp(1);exp(2);exp(2);exp(2);exp(2)];

alpha=inv(B)*C;

answer2=alpha(1)+alpha(2)*A+alpha(3)*A^2+alpha(4)*A^3+alpha(5)*A^4+alpha(6)*A^5;

P10=A;
for i=1:10
    P10=P10+A^i/factorial(i);
end

D=answer2-P10;

rho=max(abs(eig(D'*D)));

norm2=sqrt(rho);