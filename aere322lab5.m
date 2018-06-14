% Aer E 322 - Group 5 - Section 4 - Lab 5

clear,clc;

% Part (iv)

% read in strain gage data
voltage = dlmread('strain_vs_time_yellow_box.csv');

% convert to microstrain
conversion = 0.2616/1000; % V/ue
microstrain = voltage/conversion/2;

% cut off all unnecessary data
microstraincut = microstrain(5902:12000);

% Part (v)

% get values from calibration test
elongation = [0 0.00046 0.00108 0.00169 0.00231 0.00292]; % inches
load = [0.2 85.01 177.5 250.2 305 355]; % lbf

% fit a line to the data
fun = @(c,elongation) c(1)*elongation + c(2);
c = lsqcurvefit(fun, [100,0], elongation, load);
n=5021;
elongfit = linspace(0,0.003,n);
for i=1:n
    loadfit(i) = c(1)*elongfit(i) + c(2);
end

% plot the line
figure(1);
hold on;
plot(elongation,load,'*');
plot(elongfit,loadfit);
xlabel('Elongation (in)');
ylabel('Load (lbf)');
title('Load vs Elongation');

% Part (vi)

% read in the elongation data
rawdata = dlmread('Specimen_RawData_1.csv');
rawelong = rawdata(:,1);
rawload = rawdata(:,2);

% subtract off error
rawelongnew = (rawload - c(2))/c(1);
rawelongcorr = rawelong - rawelongnew;

% Part (vii)

% compute strain and area
L = 6; % inches
microstraincorr = rawelongcorr/L*1000000;
A = 0.02*2; % in^2
stress = rawload/A; % psi

% fit a line to the data
fun2 = @(c2,microstraincorr) c2(1)*microstraincorr + c2(2);
c2 = lsqcurvefit(fun2, [100,0], microstraincorr, stress);
strainfit = linspace(0,0.7,n);
for i=1:n
    stressfit(i) = c2(1)*strainfit(i) + c2(2);
end

% plot the data
figure(2);
plot(microstraincorr,stress,'*');
xlabel('Corrected Strain (10^-^6)');
ylabel('Stress (psi)');
title('Stress vs Corrected Strain Data');

% plot the line
figure(3);
plot(stressfit,strainfit*10^3);
ylabel('Strain (10^-^6)');
xlabel('Stress (psi)');
title('Stress vs Strain Fit');

% Part (viii)

% set up time vectors
for i=1:length(microstraincut);
    timesg(i) = 0.01*i;
end
for j=1:length(rawload)
    timelf(j) = 0.1*j;
end

% plot the two data sets
figure(4);
hold on;
plot(timesg,microstraincut);
plot(timelf,microstraincorr,'g');
xlabel('Time (sec)');
ylabel('Strain (10^-^6)');
title('Strain vs Time');
legend('Strain Gage','Load Frame');

% define span for data smoothing
span = 21;
reach = (span-1)/2;

% carry out first pass
newstrain = zeros(length(microstraincut)-span+1,1);
for i = 1:(length(microstraincut)-span+1)
    for j = 0:(span-1)
        newstrain(i) = newstrain(i) + microstraincut(i+j)/span;
    end
end

% clip time vector accordingly
timesg = timesg((reach+1):(length(microstraincut)-reach));

% plot the smoothed data
figure(5);
hold on;
plot(timesg,newstrain);
plot(timelf,microstraincorr,'g');
xlabel('Time (sec)');
ylabel('Strain (10^-^6)');
title('Strain vs Time');
legend('Strain Gage','Load Frame');

% clip off unnecessary flat parts of data
newstrain = newstrain(730:5750);
timesg = timesg(1:length(newstrain));

% replot the data
figure(6);
hold on;
plot(timesg,newstrain);
plot(timelf,microstraincorr,'g');
xlabel('Time (sec)');
ylabel('Strain (10^-^6)');
title('Strain vs Time');
legend('Strain Gage','Load Frame');

% Part (x)

% plot the strain vs stress
figure(7);
plot(stressfit,newstrain);
xlabel('Stress (psi)');
ylabel('Strain (10^-^6)');
title('Strain vs Stress');

% fit the data to a line
fun3 = @(c3,stressfit) c3(1)*stressfit + c3(2);
c3 = lsqcurvefit(fun3, [100,0], stressfit', newstrain);
newstrain2 = linspace(-410,-400,n);
for i=1:n
    newstrain2(i) = c3(1)*stressfit(i) + c3(2);
end

% plot the line
figure(8);
plot(stressfit,newstrain2);
xlabel('Stress (psi)');
ylabel('Strain (10^-^6)');
title('Strain vs Stress Fit');

% Part (xi)

% calculate the Young's moduli for each line
E1 = c2(1)*10^6/1000;
E2 = 1/c3(1)*10^6;