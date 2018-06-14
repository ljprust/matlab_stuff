function stiffness = stiff(E,I,x1,x2,dim,skip)

% compute length of member
L = x2 - x1;

% initialize stiffness matrix
k = zeros(dim,dim);

% fill in stiffness matrix
k(1+skip,1+skip) = 12;
k(1+skip,2+skip) = 6*L;
k(1+skip,3+skip) = -12;
k(1+skip,4+skip) = 6*L;
k(2+skip,1+skip) = 6*L;
k(2+skip,2+skip) = 4*L^2;
k(2+skip,3+skip) = -6*L;
k(2+skip,4+skip) = 2*L^2;
k(3+skip,1+skip) = -12;
k(3+skip,2+skip) = -6*L;
k(3+skip,3+skip) = 12;
k(3+skip,4+skip) = -6*L;
k(4+skip,1+skip) = 6*L;
k(4+skip,2+skip) = 2*L^2;
k(4+skip,3+skip) = -6*L;
k(4+skip,4+skip) = 4*L^2;

% if ((skip == 2)||(skip == 6))
%     E = 4E9;
%     I = 10^(-10);
% end

stiffness = k*E*I/L^3;