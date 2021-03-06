function klocal = kloc(E,I,A,L)

k = zeros(6);

k(1,1) = A*L*L/I;
k(1,4) = -k(1,1);
k(4,1) = k(1,4);
k(4,4) = k(1,1);
k(2,2) = 12;
k(2,3) = 6*L;
k(3,2) = k(2,3);
k(3,3) = 4*L*L;
k(2,5) = -12;
k(2,6) = k(2,3);
k(3,5) = -k(2,6);
k(3,6) = 2*L*L;
k(5,2) = -12;
k(5,3) = k(3,5);
k(6,2) = k(2,6);
k(6,3) = k(3,6);
k(5,5) = 12;
k(5,6) = -6*L;
k(6,5) = k(5,6);
k(6,6) = k(3,3);

klocal = E*I/L/L/L*k;