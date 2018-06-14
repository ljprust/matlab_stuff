% Logan Prust, Alexander Criswell - Phys 322L - Muon Lifetime

clear;

format long;

data=dlmread('muon_data_final.csv');

x=data(:,1);
y=data(:,2);
% percentuncert=data(:,3);

n=length(x);

% cut:
n=floor(n/2);
x=x(1:n);
y=y(1:n);

W=zeros(n,n);

for i=1:n
    X(i,1)=x(i);
    X(i,2)=1;
    if y(i)==0
        W(i,i)=0;
%     elseif percentuncert(i)>=0.25
%         W(i,i)=0;
    else
        W(i,i)=sqrt(y(i));
    end
    W2(i)=W(i,i);
end

betaW=inv(transpose(X)*W*X)*transpose(X)*W*y;

tauW=-1/betaW(1)/1000

beta=inv(transpose(X)*X)*transpose(X)*y;

tau=-1/beta(1)/1000

% beta2=lscov(x,y,W2/max(W2));

% for j=1:n
%     f(j)=betaW(2)+betaW(1)*x(j);
% end

% chi2=(transpose(y)-f)*inv(W*W)*(y-transpose(f))

% chi2_2=0;
% 
% for k=1:n
%     if y(k) ~= 0
%         chi2_2=chi2_2+(y(k)-f(k))^2/y(k);
%     end
% end
% 
% chi2_2

xbar=mean(x);
ybar=mean(y);

ssxx=0;
ssyy=0;
ssxy=0;
for i=1:n
    ssxx=ssxx+(x(i)-xbar)^2;
    ssyy=ssyy+(y(i)-ybar)^2;
    ssxy=ssxy+(x(i)-xbar)*(y(i)-ybar);
end

sW=sqrt((ssyy-betaW(1)*ssxy)/(n-2));

SEaW=sW*sqrt(1/n+xbar^2/ssxx);

SEbW=sW/sqrt(ssxx);

s=sqrt((ssyy-beta(1)*ssxy)/(n-2));

SEa=s*sqrt(1/n+xbar^2/ssxx);

SEb=s/sqrt(ssxx);

sigmatauW=SEbW/betaW(1)^2/1000

sigmatau=SEb/beta(1)^2/1000