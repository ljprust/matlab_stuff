% Logan Prust - Math 407 - Project 1 - Part B

clear,clc;
format long;

fprintf('%1s \n\n','Logan Prust - Math 407 - Project 1 - Part B');

% load L matrix from file

load('m407507proj1matrices.mat', 'L12');
L=L12;

% find dimensions of L

Lsize=size(L);
m=Lsize(1);
n=Lsize(2);

% ----------------------------------------------
% ---------------- Reduced SVD -----------------
% ----------------------------------------------

% find eigenvalues and eigenvectors of L*L

[v,lambda]=eig(L'*L);

% get rid of eigenvalues below 10^(-13)

for i=1:m
    if lambda(i,i)<10^(-13)
        lambda(i,i)=0;
    end
end

% initialize some stuff

h=1;
vcounter=ones(1,m);

% put eigenvalues in decreasing order

for k=1:m
    
    % initialize counter
    
    p=0;
    
    % initialize maximum eigenvalue
    
    maxlambda=0;
    
    % find maximum eigenvalue
    
    for i=1:m
        if lambda(i,i)>maxlambda
            maxlambda=lambda(i,i);
            j=i;
            p=1;
        end
    end
    
    % put eigenvalue in the right place
    
    lambdaordered(k)=maxlambda;
    
    % erase it from where it used to be
    
    lambda(j,j)=0;
    
    % rearrage eigenvectors accordingly
    
    if p==1
        
        vordered(:,h)=v(:,j);
        vcounter(j)=0;
        h=h+1;
        
    end
end

% find rank of L

r=h-1;

% put the rest of the eigenvectors at the end

for i=1:m
    
    if vcounter(i)==1
        
        vordered(:,h)=v(:,i);
        h=h+1;
        
    end
    
end

% find singular values

sigma=sqrt(lambdaordered);

% create Sigma hat

Sigmahat=diag(sigma(1:r));

% create V matrix

V=vordered;

for i=1:r
    
    % Vhat is just the first r columns of V
    
    Vhat(:,i)=V(:,i);
    
    % compute Uhat
    
    Uhat(:,i)=1/sigma(i)*L*V(:,i);
    
end

% output solutions to the screen

fprintf('%1s \n\n','Reduced SVD');

fprintf('%1s \n\n','Uhat =');
disp(Uhat);

fprintf('%1s \n\n','Sigma hat =');
disp(Sigmahat);

fprintf('%1s \n\n','Vhat =');
disp(Vhat);

fprintf('%1s \n\n','L - Uhat Sigma hat Vhat* =');
disp(L-Uhat*Sigmahat*Vhat');

fprintf('%1s \n\n','We can see that L - Uhat Sigma hat Vhat is effectively the zero matrix.');

% ----------------------------------------------
% ------------------ Full SVD ------------------
% ----------------------------------------------

% define U (to be modified later)

U=Uhat;

k=r;

% set a threshold for extending orthonormal basis

% (this is about the largest threshold I can have
% without having an exorbitantly large run time)

epsilon=0.8;

% implement basis extension described in Remark 4.5.3

while k<m
    for i=1:m
        
        % generate random vector
        
        w(i,1)=rand(1);
        
    end
    
    % initialize projection vector
    
    proj=zeros(m,1);
    
    % accumulate projection vector
    
    for i=1:k
        proj=proj+dot(U(:,i),w)*U(:,i);
    end
    
    % find perpendicular vector
    
    z=w-proj;
    
    % find norm of perpendicular vector
    
    squaresum=0;
    
    for j=1:m
        
        squaresum=squaresum+z(j)^2;
        
    end
    
    c=sqrt(squaresum);
    
    % only accept vector if perpendicular component is large enough
    
    if c>epsilon
        
        % normalize perpendicular vector
        
        U(:,k+1)=z/c;
        
        k=k+1;
    end
end

% end basis extension

% define Sigma (to be modified later)

Sigma=Sigmahat;

% pad Sigma with zeros

for i=r+1:m
    for j=r+1:n
        
        Sigma(i,j)=0;
        
    end
end

% output solutions to the screen

fprintf('%1s \n\n','Full SVD');

fprintf('%1s \n\n','U =');
disp(U);

fprintf('%1s \n\n','Sigma =');
disp(Sigma);

fprintf('%1s \n\n','V =');
disp(V);

fprintf('%1s \n\n','L - U Sigma V* =');
disp(L-U*Sigma*V');

fprintf('%1s \n\n','We can see that L - U Sigma V is effectively the zero matrix.');

fprintf('%1s \n\n','U*U =');
disp(U'*U);

fprintf('%1s \n\n','We can see that U*U is the identity matrix.');

fprintf('%1s \n\n','V*V =');
disp(V'*V);

fprintf('%1s \n\n','We can see that V*V is the identity matrix.');

% ----------------------------------------------
% --------------- Numerical Rank ---------------
% ----------------------------------------------

% redefine epsilon

epsilon=10^(-4);

% initialize numerical rank at zero

k0=0;

% find numerical rank using Corollary 4.3.3

for i=1:r+1
    if epsilon<sigma(i)
        
        k0=i;
        
    end
end

% output solution to the screen

fprintf('%1s \n\n','Numerical Rank k0 =');
disp(k0);

% ----------------------------------------------
% ------------ Rank-3 Approximation ------------
% ----------------------------------------------

% initialize Lk3 at zero

Lk3=0;
k=k0;

% find the best rank-3 approximation to L using Theorem 4.3.2

for i=1:k
    Lk3=Lk3+sigma(i)*U(:,i)*V(:,i)';
end

% output solution to screen

fprintf('%1s \n\n','Best Rank-3 Approximation to L =');
disp(Lk3);