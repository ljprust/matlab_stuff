clear,clc;
x = 1:10;
y = [0.84 0.81 0.87 0.67 0.58 0.71 0.67 0.68 0.68 0.47];

fun = @(c,x) c(1)+c(2).*x.^2+c(3)*x.^4+c(4)*x.^6;
c = lsqcurvefit(fun, [1,1,1,1], x, y);
% fun = @(c,x) c(1)*exp(x);
% c = lsqcurvefit(fun, [1], x, y);

c(3) = c(3)*0.6;
c(4) = 0.00000000001;

n=1100;
for i=1:n
    xfit(i) = i*0.01;
    yfit(i) = c(1)+c(2)*xfit(i)^2+c(3)*xfit(i)^4+c(4)*xfit(i)^6;
    % yfit(i) = c(1)*exp(xfit(i));
end

figure(1);
hold on;
plot(x,y,'*');
plot(xfit,yfit,'r');
axis([0 11 0 2]);
xlabel('Venue');
ylabel('Percentile');