% Logan Prust - 29 June 2016

data = dlmread('ellipsedata.csv');
alpha = data(:,1);
alphapm = data(:,2);
beta = data(:,3);
betapm = data(:,4);
f0 = data(:,5);
f0pm = data(:,6);

alphareal = -2.35771;
alphaerrreal = 0.06096713;
betareal = -0.338869;
betaerrreal = 0.21956826;
f0real = 3.37326;
f0errreal = 0.115022588;

corr_ab = -0.729;
corr_af = 0.630;
corr_bf = -0.858;

magnitude = 1;

if (magnitude == 1)
    f0 = f0*10^(-7);
    f0pm = f0pm*10^(-7);
    f0real = f0real*10^(-7);
    f0errreal = f0errreal*10^(-7);
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

[mu(1), sigma(1)] = normfit(alpha);
[mu(2), sigma(2)] = normfit(beta);
[mu(3), sigma(3)] = normfit(f0);

% sumab=0;
% sumbf=0;
% sumaf=0;
% sumaa=0;
% sumbb=0;
% sumff=0;
% 
% for j=1:n
%     sumab = sumab + (alpha(j)-mu(1))*(beta(j)-mu(2));
%     sumbf = sumbf + (beta(j)-mu(2))*(f0(j)-mu(3));
%     sumaf = sumaf + (alpha(j)-mu(1))*(f0(j)-mu(3));
%     sumaa = sumaa + (alpha(j)-mu(1))*(alpha(j)-mu(1));
%     sumbb = sumbb + (beta(j)-mu(2))*(beta(j)-mu(2));
%     sumff = sumff + (f0(j)-mu(3))*(f0(j)-mu(3));
% end
% 
Vtemp = zeros(3);
% Vtemp(1,1) = sumaa/n;
% Vtemp(2,2) = sumbb/n;
% Vtemp(3,3) = sumff/n;
% Vtemp(1,2) = sumab/n;
% Vtemp(2,1) = sumab/n;
% Vtemp(1,3) = sumaf/n;
% Vtemp(3,1) = sumaf/n;
% Vtemp(2,3) = sumbf/n;
% Vtemp(3,2) = sumbf/n;

Vtemp(1,1) = alphaerrreal*alphaerrreal;
Vtemp(2,2) = betaerrreal*betaerrreal;
Vtemp(3,3) = f0errreal*f0errreal;
Vtemp(1,2) = corr_ab*alphaerrreal*betaerrreal;
Vtemp(2,1) = Vtemp(1,2);
Vtemp(1,3) = corr_af*alphaerrreal*f0errreal;
Vtemp(3,1) = Vtemp(1,3);
Vtemp(2,3) = corr_bf*betaerrreal*f0errreal;
Vtemp(3,2) = Vtemp(2,3);

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

p = [(1-0.6827)/2 (1-0.9545)/2 (1-0.9973)/2];
quant = sqrt(2)*abs(erfinv(2*p-1));

a = sqrt(quant*lambda1);
b = sqrt(quant*lambda2);

angle = atan2(vectors(2,maxindex),vectors(1,maxindex));

t = linspace(0,2*pi,n);
x = zeros(n,3);
y = zeros(n,3);

figure(1);
hold on;
for i=1:length(p)
    X = a(i)*cos(t);
    Y = b(i)*sin(t);
    x(:,i) = param1real + X*cos(angle) - Y*sin(angle);
    y(:,i) = param2real + X*sin(angle) + Y*cos(angle);
    plot(x(:,i),y(:,i),'r-');
    % axis equal;
end

plot(param1vec,param2vec,'b*');

dataout = [param1 param2 param1real param2real n 0 0 0 0 0;param1vec param1vecpm param2vec param2vecpm x y];

dlmwrite('dataforroot.txt',dataout,'delimiter','\t','precision',6);