% Logan Prust - 20 June 2016

data = dlmread('ClusterConvHist.csv');
alpha = data(:,1);
alphapm = data(:,2);
beta = data(:,3);
betapm = data(:,4);
f0 = data(:,5);
f0 = f0*10^(-7);
f0pm = data(:,6);
f0pm = f0pm*10^(-7);

n = length(alpha);
seed = zeros(n,1);
for i=1:n
    seed(i) = i;
end

alphareal = -2.47;
betareal = -0.24;
f0real = 3.23;
f0real = f0real*10^(-7);

% figure(1);
% hold on;
% errorbar(seed,alpha,alphapm,'*');
% title('Distribution of alpha Values');
% axis([-5 205 -2.7 -2.25]);
% plot([-5 205],[-2.47 -2.47],'-');
% xlabel('Run Number (and Random Seed)');
% ylabel('alpha');
% 
% figure(2);
% hold on;
% errorbar(seed,beta,betapm,'*');
% title('Distribution of beta Values');
% axis([-5 205 -1.3 0.5]);
% plot([-5 205],[-0.24 -0.24],'-');
% xlabel('Run Number (and Random Seed)');
% ylabel('beta');
% 
% figure(3);
% hold on;
% errorbar(seed,f0,f0pm,'*');
% title('Distribution of f0 Values');
% axis([-5 205 -Inf Inf]);
% plot([-5 205],[3.23e-7 3.23e-7],'-');
% xlabel('Run Number (and Random Seed)');
% ylabel('f0');

% figure(4);
% histfit(alpha-alphareal);
% title('alpha');
% 
% figure(5);
% histfit(beta-betareal);
% title('beta');
% 
% figure(6);
% histfit(f0-f0real);
% title('f0');

[mu(1), sigma(1)] = normfit(alpha);
[mu(2), sigma(2)] = normfit(beta);
[mu(3), sigma(3)] = normfit(f0);

error(1) = 100*abs((mu(1)-alphareal)/alphareal);
error(2) = 100*abs((mu(2)-betareal)/betareal);
error(3) = 100*abs((mu(3)-f0real)/f0real);

% alphanorm = (alpha-mu(1))./alphapm;
% betanorm = (beta-mu(2))./betapm;
% f0norm = (f0-mu(3))./f0pm;
% 
% [mu2(1), sigma2(1)] = normfit(alphanorm);
% [mu2(2), sigma2(2)] = normfit(betanorm);
% [mu2(3), sigma2(3)] = normfit(f0norm);

sumab=0;
sumbf=0;
sumaf=0;
sumaa=0;
sumbb=0;
sumff=0;

for j=1:n
    sumab = sumab + (alpha(j)-mu(1))*(beta(j)-mu(2));
    sumbf = sumbf + (beta(j)-mu(2))*(f0(j)-mu(3));
    sumaf = sumaf + (alpha(j)-mu(1))*(f0(j)-mu(3));
    sumaa = sumaa + (alpha(j)-mu(1))*(alpha(j)-mu(1));
    sumbb = sumbb + (beta(j)-mu(2))*(beta(j)-mu(2));
    sumff = sumff + (f0(j)-mu(3))*(f0(j)-mu(3));
end

V = zeros(3);
V(1,1) = sumaa;
V(2,2) = sumbb;
V(3,3) = sumff;
V(1,2) = sumab;
V(2,1) = sumab;
V(1,3) = sumaf;
V(3,1) = sumaf;
V(2,3) = sumbf;
V(3,2) = sumbf;

[vectors,vals] = eig(V);

lambda = zeros(1,3);
for i=1:3
    lambda(i) = vals(i,i);
end

lambdaab = lambda(1:2);
lambdabf = lambda(2:3);
lambdaaf = [lambda(1),lambda(3)];

p = [0.159 0.023 0.002];
quant = sqrt(2)*abs(erfinv(2*p-1));

s = 1;
aab = 2*sqrt(p(1)*max(lambdaab));
abf = 2*sqrt(p(1)*max(lambdabf));
aaf = 2*sqrt(p(1)*max(lambdaaf));
bab = 2*sqrt(p(1)*min(lambdaab));
bbf = 2*sqrt(p(1)*min(lambdabf));
baf = 2*sqrt(p(1)*min(lambdaaf));

if (lambda(1) > lambda(2))
    oneortwo = 1;
else
    oneortwo = 2;
end

if (lambda(2) > lambda(3))
    twoorthree = 2;
else
    twoorthree = 3;
end

if (lambda(1) > lambda(3))
    oneorthree = 1;
else
    oneorthree = 3;
end

angleab = atan2(vectors(2,oneortwo),vectors(1,oneortwo));
anglebf = atan2(vectors(3,twoorthree),vectors(2,twoorthree));
angleaf = atan2(vectors(3,oneorthree),vectors(1,oneorthree));

figure(7);
hold on;
t = linspace(0,2*pi);
X = aaf*cos(t);
Y = baf*sin(t);
w = angleaf;
x = mu(1) + X*cos(w) - Y*sin(w);
y = mu(3) + X*sin(w) + Y*cos(w);
plot(x,y,'r-');
plot(alpha,f0,'b*');
% axis equal;

