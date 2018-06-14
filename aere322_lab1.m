% Aer E 322 - Lab 1 - Section 4 - Group 5

clear,clc;

% read in data from file
data = dlmread('Lab_1_Part_3_Group_5-2.txt');

% separate data into relevant components
time = data(:,1);
pos1 = data(:,2);
pos2 = data(:,6);
pos3 = data(:,10);

[n,m] = size(time);

% input displacement sensor data
needlepos = [-.9398 1.2192 1.9939 1.4986 -.6096 -1.9177 -1.7399 -.9398 .0381 1.9431 1.8669 ...
    -.1524 -1.5875 -2.0193 -1.4732 0.3937 1.4224 2.0955 1.1049 -1.1049 -2.0066 -1.4478 ...
    -.4699 .6096 1.9177 1.5748 -.4572 -1.651 -1.8161 -1.143 -.2286 1.7399 0.0 .7493 -1.3589 ...
    -1.9939 -1.0033 -.9525 1.0033 1.9685 1.2319 -.8636 -1.905 -1.6256 -.7239 .3175 1.8796 ...
    1.6129 -.3175 -1.6637 -1.905 -1.27 -.2413 1.8034 1.9177 .7493 -1.3843 -2.0701 -1.0287 ...
    .127 1.1811 2.0828 1.27 -.9271 -2.0574 -1.6002 -.5842 .4318 1.9177 1.651 -.4064 -1.778 ...
    -1.9685 -1.2573 -.2413];
[junk,needlesize] = size(needlepos);

% set up time points for displacement sensor
needletime = zeros(1,needlesize);
for i = 1:needlesize
    needletime(i) = i*0.2;
end

% set span for first pass
span = 9;
reach = (span-1)/2;

% initialize new variables
newpos1 = zeros(n-span+1,1);
newpos2 = zeros(n-span+1,1);
newpos3 = zeros(n-span+1,1);

% run first pass
for i = 1:(n-span+1)
    for j = 0:(span-1)
        newpos1(i) = newpos1(i) + pos1(i+j)/span;
        newpos2(i) = newpos2(i) + pos2(i+j)/span;
        newpos3(i) = newpos3(i) + pos3(i+j)/span;
    end
end

% set span for second pass
span2 = 9;
reach2 = (span2-1)/2;

% initialize new variables
newpos1_2 = zeros(n-span-span2+2,1);
newpos2_2 = zeros(n-span-span2+2,1);
newpos3_2 = zeros(n-span-span2+2,1);

% run second pass
for i = 1:(n-span-span2+2)
    for j = 0:(span2-1)
        newpos1_2(i) = newpos1_2(i) + newpos1(i+j)/span2;
        newpos2_2(i) = newpos2_2(i) + newpos2(i+j)/span2;
        newpos3_2(i) = newpos3_2(i) + newpos3(i+j)/span2;
    end
end

% cut time data to correct length
newtime = time((reach+1):(n-reach));
newtime2 = time((reach+reach2+1):(n-reach-reach2));

% set up curve fit options
options = optimset('TolFun',1e-10);

% define fit function
fun = @(x,newtime2) x(1)*sin(x(2)*newtime2 + x(3)) + x(4);

% carry out fits
x1 = lsqcurvefit(fun, [1,2*2*pi,pi,4], newtime2, 100000*newpos1_2);
x2 = lsqcurvefit(fun, [1,8/5*2*pi,0,-2.5], newtime2, 1000*newpos2_2);
x3 = lsqcurvefit(fun, [1,1.4*2*pi,0,6.5], newtime2, 1000*newpos3_2);

% initialize fit curves
sinwave1 = zeros(1,length(newtime2));
sinwave2 = zeros(1,length(newtime2));
sinwave3 = zeros(1,length(newtime2));

% calculate fit curves
for k = 1:length(newtime2)
    sinwave1(k) = (x1(1)*sin(x1(2)*newtime2(k) + x1(3)) + x1(4))/100000;
    sinwave2(k) = (x2(1)*sin(x2(2)*newtime2(k) + x2(3)) + x2(4))/1000;
    sinwave3(k) = (x3(1)*sin(x3(2)*newtime2(k) + x3(3)) + x3(4))/1000;
end

% blahtime = newtime2(500:1000);
% blahpos = 1000*newpos2_2(500:1000);
% blahfun = @(x,blahtime) x(1)*sin(x(2)*blahtime + x(3)) + x(4);
% x = lsqcurvefit(blahfun, [1,8/5*2*pi,0,-2.5], blahtime, blahpos);
% for i = 1:length(blahtime)
%     blahwave(i) = x(1)*sin(x(2)*blahtime(i) + x(3)) + x(4);
% end
% 
% figure(100);
% hold on;
% plot(blahtime,blahpos,'b');
% plot(blahtime,blahwave,'g');

% create all relevant plots

% figure(1);
% plot(time,pos1);
% title('Run 1 Raw Data (Motion Sensor)');
% xlabel('Time (sec)');
% ylabel('Position (m)');
% 
% figure(2);
% plot(time,pos2);
% title('Run 2 Raw Data');
% xlabel('Time (sec)');
% ylabel('Position (m)');
% 
% figure(3);
% plot(time,pos3);
% title('Run 3 Raw Data');
% xlabel('Time (sec)');
% ylabel('Position (m)');
% 
% figure(4);
% plot(needletime,needlepos);
% title('Run 1 Raw Data (Needle)');
% xlabel('Time (sec)');
% ylabel('Position (m)');

figure(5);
plot(newtime2,newpos1_2);
title('Run 1 Smoothed Data');
xlabel('Time (sec)');
ylabel('Position (m)');

% figure(6);
% plot(newtime2,newpos2_2);
% title('Run 2 Smoothed Data');
% xlabel('Time (sec)');
% ylabel('Position (m)');
% 
% figure(7);
% plot(newtime2,newpos3_2);
% title('Run 3 Smoothed Data');
% xlabel('Time (sec)');
% ylabel('Position (m)');
% 
% figure(8);
% hold on;
% plot(newtime2,newpos1_2,'b');
% plot(newtime2,sinwave1,'g');
% title('Run 1 Data With Fitted Curve');
% xlabel('Time (sec)');
% ylabel('Position (m)');
% 
% figure(9);
% hold on;
% plot(newtime2,newpos2_2,'b');
% plot(newtime2,sinwave2,'g');
% title('Run 2 Data With Fitted Curve');
% xlabel('Time (sec)');
% ylabel('Position (m)');
% 
% figure(10);
% hold on;
% plot(newtime2,newpos3_2,'b');
% plot(newtime2,sinwave3,'g');
% title('Run 3 Data With Fitted Curve');
% xlabel('Time (sec)');
% ylabel('Position (m)');