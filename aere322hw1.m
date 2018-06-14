% Logan Prust - Aer E 322 - Section 4 - HW 1

clear,clc;

% read data from file
data = dlmread('HW1_Data.txt');

% separate data
time = data(:,1); % microseconds
amp = data(:,2);

% plot raw data
figure(1);
plot(time,amp);
xlabel('Time (microseconds)');
ylabel('Amplitude');
title('Raw Data');

% extract number of data points
[n,m] = size(data);

% set span for first pass
span = 11;
reach = (span-1)/2;

% carry out first pass
newamp = zeros(n-span+1,1);
for i = 1:(n-span+1)
    for j = 0:(span-1)
        newamp(i) = newamp(i) + amp(i+j)/span;
    end
end

% set span for second pass
span2 = 35;
reach2 = (span2-1)/2;

% carry out second pass
newamp2 = zeros(n-span-span2+2,1);
for i = 1:(n-span-span2+2)
    for j = 0:(span2-1)
        newamp2(i) = newamp2(i) + newamp(i+j)/span2;
    end
end

% clip edges of time data
newtime = time((reach+1):(n-reach));
newtime2 = time((reach+reach2+1):(n-reach-reach2));

% plot smoothed data
figure(2);
plot(newtime2,newamp2);
xlabel('Time (microseconds)');
ylabel('Amplitude');
title('Smoothed Data');

% define fit function
fun = @(x,newtime2)x(1)*sin(x(2)*newtime2 + x(3));

% carry out fit
x = lsqcurvefit(fun, [150,0.05*2*pi,0], newtime2, newamp2);

% create fitted curve
for k = 1:length(newtime2)
    sinwave(k) = x(1)*sin(x(2)*newtime2(k) + x(3));
end

% plot fitted curve
figure(3);
plot(newtime2,sinwave);
xlabel('Time (microseconds)');
ylabel('Amplitude');
title('Fitted Curve (sine)');

% overlay plots
figure(4);
hold on;
plot(time,amp,'b');
plot(newtime2,newamp2,'g');
plot(newtime2,sinwave,'r');
xlabel('Time (microseconds)');
ylabel('Amplitude');
title('Raw Data, Smoothed Data, and Fitted Curve');

% plot final signal
figure(6);
plot(newtime2,abs(newamp2-sinwave'));
axis([8 16 0 50]);
xlabel('Time (microseconds)');
ylabel('Amplitude');
title('Final Signal');

% define cosine fit function
fun2 = @(x2,newtime2)x2(1)*cos(x2(2)*newtime2 + x2(3) - pi/2);
% carry out cosine fit
x2 = lsqcurvefit(fun2, [150,0.05*2*pi,pi/2], newtime2, newamp2);

% create cosine fitted function
for k = 1:length(newtime2)
    coswave(k) = x2(1)*cos(x2(2)*newtime2(k) + x2(3) - pi/2);
end

% plot cosine fitted function
figure(7);
plot(newtime2,coswave);
xlabel('Time (microseconds)');
ylabel('Amplitude');
title('Fitted Curve (cosine)');