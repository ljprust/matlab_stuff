% Logan Prust - Math 407 - Project 1 - Part A

clear,clc;
format long;

fprintf('%1s \n\n','Logan Prust - Math 407 - Project 1 - Part A');

% load H matrix from file

load('m407507proj1matrices.mat', 'H12');
H=H12;

% find dimensions of H

Hsize=size(H);
m=Hsize(1);
n=Hsize(2);

% ----------------------------------------------
% ---------- Reduced QR Decomposition ----------
% ----------------------------------------------

% define W for use in Algorithm 2

W=H;

% implement Algorithm 2 from HLA

for i=1:n
    
    % find 2-norm of w vector
    
    squaresum=0;
    
    for j=1:m
        squaresum=squaresum+W(j,i)^2;
    end
    
    r(i,i)=sqrt(squaresum);
    
    % construct u vector
    
    U(:,i)=W(:,i)/r(i,i);
    
    for j=i+1:n
        
        % construct r matrix
        
        r(i,j)=dot(W(:,j),U(:,i));
        
        % redefine W matrix
        
        W(:,j)=W(:,j)-r(i,j)*U(:,i);
    end
end

% end Algorithm 2

% redefine U and r as Qhat and Rhat

Qhat=U;
Rhat=r;

% output solution to screen

fprintf('%1s \n\n','Reduced QR Decomposition');

fprintf('%1s \n\n','Qhat =');
disp(Qhat);

fprintf('%1s \n\n','Rhat =');
disp(Rhat);

fprintf('%1s \n\n','H - Qhat Rhat =');
disp(H-Qhat*Rhat);

fprintf('%1s \n\n','We can see that H - Qhat Rhat is effectively the zero matrix.');

% ----------------------------------------------
% ----------- Full QR Decomposition ------------
% ----------------------------------------------

% define Q and R (to be modified later)

Q=Qhat;
R=Rhat;

k=n;

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
        proj=proj+dot(Q(:,i),w)*Q(:,i);
    end
    
    % find perpendicular vector
    
    v=w-proj;
    
    % find norm of perpendicular vector
    
    squaresum=0;
    
    for j=1:m
        squaresum=squaresum+v(j)^2;
    end
    
    c=sqrt(squaresum);
    
    % only accept vector if perpendicular component is large enough
    
    if c>epsilon
        
        % normalize perpendicular vector
        
        Q(:,k+1)=v/c;
        
        % pad R matrix with zeros
        
        for i=1:n
            R(k+1,i)=0;
        end
        
        k=k+1;
    end
end

% end basis extension

% output solution to screen

fprintf('%1s \n\n','Full QR Decomposition');

fprintf('%1s \n\n','Q =');
disp(Q);

fprintf('%1s \n\n','R =');
disp(R);

fprintf('%1s \n\n','H - QR =');
disp(H-Q*R);

fprintf('%1s \n\n','We can see that H - QR is effectively the zero matrix.');

fprintf('%1s \n\n','Q*Q =');
disp(Q'*Q);

fprintf('%1s \n\n','We can see that Q*Q is the identity matrix.');

% ----------------------------------------------
% ----------- Least-Squares Solution -----------
% ----------------------------------------------

% input b vector

bfixed=[1;-2;0;3;1;1];

% redefine matrices for use in backsolver

A=Rhat;
b=Qhat'*bfixed;

% create augmented matrix

C=[A b];

% normalize so that diagonal elements are 1

for i=1:n
    C(i,:)=C(i,:)/C(i,i);
end

% solve for all x and accumulate in b

for i=n-1:-1:1
    for j=i+1:n
        C(i,n+1)=C(i,n+1)-C(i,j)*C(j,n+1);
    end
end

% separate solution from rest of matrix

x0=C(:,n+1);

% use Matlab backsolver to check solution

x0check=H\bfixed;

% output solution to screen

fprintf('%1s \n\n','Least Squares Solution');

fprintf('%1s \n\n','x0 =');
disp(x0);