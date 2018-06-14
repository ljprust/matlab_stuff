function stiffness2 = stiff2(E,I,A,L,t,dim,b,e)

start1 = 3*b-3;
start2 = 3*e-3;

% initialize stiffness matrix
k = zeros(dim,dim);

% fill in stiffness matrix
k(1+start1,1+start1) = A*L*L/I*cos(t)^2+12*sin(t)^2;
k(1+start1,2+start1) = (A*L*L/I-12)*cos(t)*sin(t);
k(1+start1,3+start1) = -6*L*sin(t);
k(1+start1,1+start2) = -k(1+start1,1+start1);
k(1+start1,2+start2) = -k(1+start1,2+start1);
k(1+start1,3+start2) = k(1+start1,3+start1);
k(2+start1,1+start1) = k(1+start1,2+start1);
k(2+start1,2+start1) = A*L*L/I*sin(t)^2+12*cos(t)^2;
k(2+start1,3+start1) = 6*L*cos(t);
k(2+start1,1+start2) = -k(2+start1,1+start1);
k(2+start1,2+start2) = -k(2+start1,2+start1);
k(2+start1,3+start2) = k(2+start1,3+start1);
k(3+start1,1+start1) = k(1+start1,3+start1);
k(3+start1,2+start1) = k(2+start1,3+start1);
k(3+start1,3+start1) = 4*L*L;
k(3+start1,1+start2) = -k(1+start1,3+start1);
k(3+start1,2+start2) = -k(3+start1,2+start1);
k(3+start1,3+start2) = 2*L*L;
k(1+start2,1+start1) = k(1+start1,1+start2);
k(1+start2,2+start1) = k(2+start1,1+start2);
k(1+start2,3+start1) = k(3+start1,1+start2);
k(1+start2,1+start2) = k(1+start1,1+start1);
k(1+start2,2+start2) = k(1+start1,2+start1);
k(1+start2,3+start2) = -k(1+start1,3+start2);
k(2+start2,1+start1) = k(1+start1,2+start2);
k(2+start2,2+start1) = k(2+start1,2+start2);
k(2+start2,3+start1) = k(3+start1,2+start2);
k(2+start2,1+start2) = k(1+start2,2+start2);
k(2+start2,2+start2) = k(2+start1,2+start1);
k(2+start2,3+start2) = k(3+start1,2+start2);
k(3+start2,1+start1) = k(1+start1,3+start2);
k(3+start2,2+start1) = k(2+start1,3+start2);
k(3+start2,3+start1) = k(3+start1,3+start2);
k(3+start2,1+start2) = k(1+start2,3+start2);
k(3+start2,2+start2) = k(2+start2,3+start2);
k(3+start2,3+start2) = k(3+start1,3+start1);

stiffness2 = k*E*I/L^3;