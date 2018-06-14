% Logan Prust - Aer E 322 - Section 4 - Homework 3

clear,clc;

% read in truss parameters from input file
data = dlmread('aere322hw3.txt');

% sort data
npins = data(1);
nmem = data(2);
E = data(3);
I = data(4);
A = data(5);
x = data(6:(5+npins));
y = data((6+npins):(5+2*npins));
freeflag = data((6+2*npins):(5+5*npins));

% set up matrix with free degrees of freedom
nfree=0;
for m=1:npins*3
    if (freeflag(m) == 1)
        nfree = nfree+1;
        freemat(nfree) = m;
    end
end

% sort more of the data
Pfree = data((6+5*npins):(5+5*npins+nfree));
con = data((6+5*npins+nfree):(5+5*npins+nfree+2*nmem));

% get connectivity information
for i=1:nmem
    b(i) = con(2*i-1);
    e(i) = con(2*i);
end

% initialize some stuff
dim = 3*npins;
ktotal = zeros(dim,dim);
k = zeros(dim,dim,nmem);
T = zeros(6,6,nmem);
klocal = zeros(6,6,nmem);

% calculate member properties and stiffness matrices
for i=1:nmem
    L(i) = sqrt((x(b(i))-x(e(i)))^2+(y(b(i))-y(e(i)))^2);
    t(i) = asin((y(e(i))-y(b(i)))/L(i));
    klocal(:,:,i) = kloc(E,I,A,L(i));
    T(:,:,i) = trans(t(i));
    k(:,:,i) = stiff2(E,I,A,L(i),t(i),dim,b(i),e(i));
    ktotal = ktotal + k(:,:,i);
end

% find final stiffness matrix (with free DOFs)
kfin = zeros(nfree);
for i=1:nfree
    for j=1:nfree
        kfin(i,j) = ktotal(freemat(i),freemat(j));
    end
end

% calculate displacements
dfree = kfin\Pfree;

% sort displacements by member
for i=1:nmem
    for m=1:3
        for j=1:nfree
            if freemat(j) == (3*b(i)-3+m)
                v(m,i) = dfree(j);
                break;
            else
                v(m,i) = 0;
            end
        end
    end
    for m=1:3
        for j=1:nfree
            if freemat(j) == (3*e(i)-3+m)
                v(m+3,i) = dfree(j);
                break;
            else
                v(m+3,i) = 0;
            end
        end
    end
end

% calculate the local displacements, local loads, global loads, axial 
% stresses, and axial strains
for i=1:nmem
    u(:,i) = T(:,:,i)*v(:,i);
    Q(:,i) = klocal(:,:,i)*u(:,i);
    sigma_a(i) = Q(4,i)/A;
    eps_a(i) = sigma_a(i)/E;
    F(:,i) = T(:,:,i)'*Q(:,i);
end

% Lab 9 Portion

% define c and gage positions on members
c = 0.125;
gaugepos = [0 3.5 4 3.2 0 4.25 4.75 3.5 9 7 7.25 7.25];

% calculate bending stresses and total stresses and strains
for i=1:nmem
    sigma_b(i) = -c/I*(gaugepos(i)*Q(2,i)-Q(3,i));
    sigmatot(i) = sigma_a(i)+sigma_b(i);
    epstot(i) = sigmatot(i)/E;
end

% HW 3 Portion

% initialize some stuff
iter = 1000;
c=0.125;
xvector = linspace(0,L(5),iter);

% calculate stresses and strains on surfaces of members 5 and 6
for i=1:iter
    sigma5top(i) = -c/I*(xvector(i)*Q(2,5)-Q(3,5))+sigma_a(5);
    sigma6top(i) = -c/I*(xvector(i)*Q(2,6)-Q(3,6))+sigma_a(6);
    sigma5bottom(i) = c/I*(xvector(i)*Q(2,5)-Q(3,5))+sigma_a(5);
    sigma6bottom(i) = c/I*(xvector(i)*Q(2,6)-Q(3,6))+sigma_a(6);
    eps5top(i) = sigma5top(i)/E;
    eps5bottom(i) = sigma5bottom(i)/E;
    eps6top(i) = sigma6top(i)/E;
    eps6bottom(i) = sigma6bottom(i)/E;
end

% plot results
figure(1);
hold on;
plot(xvector,eps5top*10^6);
plot(xvector,eps5bottom*10^6,'r');
legend('Top','Bottom');
xlabel('Location Along Member (in)');
ylabel('Total Strain (microns)');
title('Strain on Member 5');

figure(2);
hold on;
plot(xvector,eps6top*10^6);
plot(xvector,eps6bottom*10^6,'r');
legend('Top','Bottom');
xlabel('Location Along Member (in)');
ylabel('Total Strain (microns)');
title('Strain on Member 6');