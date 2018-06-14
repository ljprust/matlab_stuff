% Logan Prust - Aer E 322 - Homework 2

clear,clc;

% read in the beam parameters from the input file
data = dlmread('aere322hw2redo.txt');
nmem = data(1);
npins = nmem+1;
E = data(2);
I = data(3);

% read in more stuff
i = 1;
j = 1;
m = 1;
while i<=(npins*7)
    % like the node locations, BCs, and applied loads
    x(m) = data(3+i);
    u(j) = data(4+i); % vertical displacements
    freeflag(j) = data(5+i);
    u(j+1) = data(6+i); % rotational displacements
    freeflag(j+1) = data(7+i);
    F(j) = data(8+i); % applied forces
    F(j+1) = data(9+i); % applied moments
    i=i+7;
    j=j+2;
    m=m+1;
end

% construct the stiffness matrices
dim = 2*nmem+2;
k = zeros(dim,dim,nmem);
ktotal = zeros(dim);
for j=1:nmem
    % find local matrix
    k(:,:,j) = stiff(E,I,x(j),x(j+1),dim,j*2-2);
    % add them together
    ktotal = ktotal + k(:,:,j);
end

% set up matrix with free degrees of freedom
nfree=0;
for m=1:npins*2
    if (freeflag(m) == 1)
        nfree = nfree+1;
        freemat(nfree) = m;
    end
end

% set up matrix with applied loads on free DOFs
P = zeros(nfree,1);
for i=1:nfree
    P(i) = F(freemat(i));
end

% find final stiffness matrix (with free DOFs)
kfin = zeros(nfree);
for i=1:nfree
    for j=1:nfree
        kfin(i,j) = ktotal(freemat(i),freemat(j));
    end
end

% calculate free DOFs
ufree = inv(kfin)*P;

% integrate free DOFs with givens
for i=1:nfree
    u(freemat(i)) = ufree(i);
end

% find reaction loads
Q = ktotal*u';