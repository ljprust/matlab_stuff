% Logan Prust - 29 June 2016

data = dlmread('ellipsedata.csv');
alpha = data(:,1);
alphapm = data(:,2);
beta = data(:,3);
betapm = data(:,4);
f0 = data(:,5);
f0pm = data(:,6);

alphareal = -2.35771;
betareal = -0.338869;
f0real = 3.37326;
alphareal0 = -2.47;
betareal0 = -0.24;
f0real0 = 3.23;

magnitude = 1;

if (magnitude == 1)
    f0 = f0*10^(-7);
    f0pm = f0pm*10^(-7);
    f0real = f0real*10^(-7);
    f0real0 = f0real0*10^(-7);
end

param1 = input('Parameter 1: ');
param2 = input('Parameter 2: ');

if (param1 == 1)
    param1vec = alpha;
    param1vecpm = alphapm;
    param1real = alphareal;
elseif (param1 == 2)
    param1vec = beta;
    param1vecpm = betapm;
    param1real = betareal;
elseif (param1 == 3)
    param1vec = f0;
    param1vecpm = f0pm;
    param1real = f0real;
end

if (param2 == 1)
    param2vec = alpha;
    param2vecpm = alphapm;
    param2real = alphareal;
elseif (param2 == 2)
    param2vec = beta;
    param2vecpm = betapm;
    param2real = betareal;
elseif (param2 == 3)
    param2vec = f0;
    param2vecpm = f0pm;
    param2real = f0real;
end

n = length(alpha);

% [mu(1), sigma(1)] = normfit(alpha);
% [mu(2), sigma(2)] = normfit(beta);
% [mu(3), sigma(3)] = normfit(f0);

sumab = (alphareal-alphareal0)*(betareal-betareal0);
sumbf = (betareal-betareal0)*(f0real-f0real0);
sumaf = (alphareal-alphareal0)*(f0real-f0real0);
sumaa = (alphareal-alphareal0)*(alphareal-alphareal0);
sumbb = (betareal-betareal0)*(betareal-betareal0);
sumff = (f0real-f0real0)*(f0real-f0real0);

Vtemp = zeros(3);
Vtemp(1,1) = sumaa;
Vtemp(2,2) = sumbb;
Vtemp(3,3) = sumff;
Vtemp(1,2) = sumab;
Vtemp(2,1) = sumab;
Vtemp(1,3) = sumaf;
Vtemp(3,1) = sumaf;
Vtemp(2,3) = sumbf;
Vtemp(3,2) = sumbf;

V = zeros(2);
V(1,1) = Vtemp(param1,param1);
V(2,2) = Vtemp(param2,param2);
V(1,2) = Vtemp(param1,param2);
V(2,1) = V(1,2);

[vectors,vals] = eig(V);

lambda = zeros(1,2);
for i=1:2
    lambda(i) = vals(i,i);
end

lambda1 = max(lambda);
lambda2 = min(lambda);

if (lambda(1) > lambda(2))
    maxindex = 1;
else
    maxindex = 2;
end

p = [0.159 0.023 0.002];
quant = sqrt(2)*abs(erfinv(2*p-1));

a = 2*sqrt(quant*lambda1);
b = 2*sqrt(quant*lambda2);

angle = atan2(vectors(2,maxindex),vectors(1,maxindex));

t = linspace(0,2*pi,n);
x = zeros(n,3);
y = zeros(n,3);

figure(1);
hold on;
for i=1:length(p)
    X = a(i)*cos(t);
    Y = b(i)*sin(t);
    x(:,i) = mu(param1) + X*cos(angle) - Y*sin(angle);
    y(:,i) = mu(param2) + X*sin(angle) + Y*cos(angle);
    plot(x(:,i),y(:,i),'r-');
    % axis equal;
end

plot(param1vec,param2vec,'b*');

% dataout = [param1 param2 param1real param2real n 0;param1vec param1vecpm param2vec param2vecpm x' y'];

% dlmwrite('dataforroot.txt',dataout,'delimiter','\t','precision',6);